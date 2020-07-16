import sys
import pyranges as pr
from cyvcf2 import VCF, Writer
import time
import gzip
import pandas as pd
import argparse
import multiprocessing
import io
from collections import defaultdict


########################################
##          Global Variables          ##
########################################


extra_annotation_choices = ["all", "mf", "best", "pops", "AFR", "AMR", "EAS", "EUR", "OTH", "SAS", "full", "mis"]


########################################
##          Parser Arguments          ##
########################################

def add_annotation(parser):

    p = parser.add_parser("annotate",
        help="Annotate SV VCF File",
        description="Annotate your SV VCF files with Population Allele Frequency Info"
    )

    req = p.add_argument_group("Required Arguments")
    req_one = p.add_argument_group("Requires Only One Argument")
    opt = p.add_argument_group("Optional Arguments")
    p._optionals.title = "help"

    req.add_argument("-v", "--vcf",
                     metavar="INPUT VCF",
                     required = True,
                     help="Path and/or name of the VCF file to annotate")

    req.add_argument("-o", "--out",
                     metavar="OUTPUT VCF",
                     required = True,
                     help="Path and/or name of the output VCF file")

    req_one.add_argument("-b", "--bed",
                     metavar="SOURCE BED",
                     help="Path and/or name of the combined sources AF bed file. (--pickled-source can be used instead of --bed) (NOTE: --bed takes high priority then --pickled-source")

    req_one.add_argument("-p", "--pickled-source",
                         metavar = "Pickled Source Data",
                         help = "Path and/or name of the pickled source data to use. (--bed can be used instead of --pickled-source) (NOTE: if --bed is provided, --pickled-source is ignored.)")

    opt.add_argument("-f", "--minf",
                   metavar="MINIMUM OVERLAP FRACTION",
                   nargs = "*",
                   help=("A space sperated list of minimum reciprocal overlap fractions required between SVs for each source"
                         " listed with the `-s` option. If `-s` is not used, only the first minf will be used and will be applied"
                         " to all sources. minf values must be between 0.0 and 1.0. (Defualt = 0.001)"))

    opt.add_argument("-s", "--sources",
                   metavar="SOURCES TO ANNOTATE",
                   nargs = "*",
                   help=("Space sperated list of data sources to use for annotation. If '-s' is not used, all sources available in the"
                         " source bed file will be used. (Example: ' -s CCDG CEPH gnomAD 1000G_Smoove ' )"))

    opt.add_argument("-a", "--ann",
                   metavar="EXTRA ANNOTATIONS",
                   nargs = "*",
                   choices = extra_annotation_choices,
                   help=("By default, only the maxAF, max Hets and Homalt counts, and max PopAF are annotated in the output VCF file."
                         " `-a` can be used to add additional anntations, with each anntotation seperated by a space."
                         " (Example ' -a mf best pops ' ). Choices = [{}]").format(", ".join(extra_annotation_choices)))

    opt.add_argument("-ci", "--ci",
                   metavar="USE CI BOUNDARIES",
                   choices=["in","out"],
                   help="If argument selected, use 'inner' or 'outer' confidence intervals (CIPOS95, CIEND95) for SV boundaries. Choices = [in, out]")

    opt.add_argument("--cpu",
                   metavar="CPU Count",
                   default = 1,
                   help="The number of cpus to use for multi-threading (Default = 1)")

    p.set_defaults(func=annotate)

########################################
##              Functions             ##
########################################

def process_bed_source(bed_file):

    ## process bed file containing SV AFs
    ## save info into datas in source specific manner
    ## create source specific lists for use in pyranges

    sources = set()
    datas = defaultdict(lambda: defaultdict(list))
    bed_lists = defaultdict(lambda: defaultdict(list))
    bed_headers = []

    bed = gzip.open(bed_file, "rt", encoding="utf-8") if bed_file.endswith(".gz") else io.open(bed_file, "rt", encoding="utf-8")
    print("\nReading the following data source file: {}".format(bed_file))

    header = []
    header_found = False
    features = ["chrom","start","end","svtype","sv_id"]
    for line in bed:

        ## get header
        if line.startswith("#"):
            header = line.strip().replace("#","").split("\t")
            header_found = True
            bed_headers.append(header[3])
            bed_headers.extend(header[6:len(header)])
            continue

        ## Check for header
        if header_found == False:
            raise ValueError("A header is required in the bed source file")

        ## Set up data
        line_dict = dict(zip(header, line.strip().split("\t")))
        source = line_dict["SOURCE"]
        sv_id = line_dict["SV_ID"]

        ## Collect info
        datas[source][sv_id].append(line_dict["SVTYPE"])
        datas[source][sv_id].extend([line_dict[field] for field in header[6:len(header)]])
    
        ## add freature info for
        for i,key in enumerate(features):
            bed_lists[source][features[i]].append(line_dict[key.upper()])

        ## add sources sources 
        sources.add(source)

    bed.close()

    return(souces,datas,bed_lists, bed_headers)


def process_pickled_source(pickled_source):

    import pickle

    print("\nReading pickled data source")
    pickle_file = open(pickled_source,"rb")
    pickled_data = pickle.load(pickle_file)
    pickle_file.close()
    
    return(pickled_data["sources"], pickled_data["datas"], pickled_data["bed_lists"], pickled_data["bed_headers"])


def get_header_types(bed_headers):

    header_cols = dict(zip(bed_headers, range(len(bed_headers))))
    header_types = {}
    for col in bed_headers:
        if col not in header_types:
            header_types[col] = "NA"
        if col == "AF":
            field_type = "Float"
            header_types[col] = field_type
        elif "_AF" in col:
            field_type = "Float"
            header_types[col] = field_type
        else:
            field_type = "Integer"
            header_types[col] = field_type

    return(header_types,header_cols)


def get_overlap(groupby_list):
    ## find the amount of overlap between two sets of coordinates
    ## expects lists comprising start and end coordinates
    ## returns list of overlaps for each set of coordinates
    ## used by convert_dict function
    chunk_list = []
    for pos_list in groupby_list:
        s1 = pos_list[0]
        e1 = pos_list[1]
        s2 = pos_list[2]
        e2 = pos_list[3]
        size = int(e1) - int(s1)
        chunk = 0
        if int(s1) <= int(s2):
            if int(e1) <= int(e2):
                chunk = int(e1) - int(s2)
            elif int(e1) > int(e2):
                chunk = int(e2) - int(s2)
        elif int(s1) > int(s2):
            if int(e1) <= int(e2):
                chunk = int(e1) - int(s1)
            elif int(e1) > int(e2):
                chunk = int(e2) - int(s1)
        chunk_list.append(chunk)
    return(chunk_list)


def get_fraction1(groupby_list):
    ## similar to get_overlap function with same logic and requirements
    ## finds overlap fraction between overlap and first set of coordinates
    ## returns list of overlap fractions for each set of coordinates
    fraction_list = []
    for pos_list in groupby_list:
        s1 = pos_list[0]
        e1 = pos_list[1]
        s2 = pos_list[2]
        e2 = pos_list[3]
        size = int(e1) - int(s1)
        chunk = 0
        if int(s1) <= int(s2):
            if int(e1) <= int(e2):
                chunk = int(e1) - int(s2)
            elif int(e1) > int(e2):
                chunk = int(e2) - int(s2)
        elif int(s1) > int(s2):
            if int(e1) <= int(e2):
                chunk = int(e1) - int(s1)
            elif int(e1) > int(e2):
                chunk = int(e2) - int(s1)
        fraction_list.append(float(chunk)/size)
    return(fraction_list)


def get_fraction2(groupby_list):
    ## same as get_fraction1 function
    ## finds overlap fraction between overlap and second set of coordinates
    ## returns list of overlap fractions for each set of coordinates
    fraction_list = []
    for pos_list in groupby_list:
        s1 = pos_list[2]
        e1 = pos_list[3]
        s2 = pos_list[0]
        e2 = pos_list[1]
        size = int(e1) - int(s1)
        chunk = 0
        if int(s1) <= int(s2):
            if int(e1) <= int(e2):
                chunk = int(e1) - int(s2)
            elif int(e1) > int(e2):
                chunk = int(e2) - int(s2)
        elif int(s1) > int(s2):
            if int(e1) <= int(e2):
                chunk = int(e1) - int(s1)
            elif int(e1) > int(e2):
                chunk = int(e2) - int(s1)
        fraction_list.append(float(chunk)/size)
    return(fraction_list)


def convert_dict(pr):
    ## change a pyranges object into a dictionary
    ## input is a pyranges object
    ## returns a dictionary with main key, "SV_ID"
    df = pr.as_df()
    df["Overlap"] = df[["Start","End","Start_b","End_b","SV_ID_b"]].values.tolist()
    df["Fraction"] = df["Overlap"]
    df["Fraction_b"] = df["Overlap"]
    my_dict = df.groupby("SV_ID").agg({"Start": list, \
                                       "End": list, \
                                       "Start_b": list, \
                                       "End_b": list, \
                                       "SV_ID_b": list, \
                                       "Overlap": get_overlap, \
                                       "Fraction": get_fraction1, \
                                       "Fraction_b": get_fraction2 }).transpose().to_dict()
    return(my_dict)


def reciprocal_overlap(my_dict,source,minfs):
    ## using minf, filter products of pyranges.join based on amount of reciprocal overlap
    ## expect dictionary (from convert_dict)
    ## returns dictionary of those that pass the filter

    minf = float(minfs[source])
    pass_filter = defaultdict(lambda: defaultdict(list))

    for sv_id in my_dict:

        for i in range(len(my_dict[sv_id]["Fraction"])):
            start = int(my_dict[sv_id]["Start"][i])
            end = int(my_dict[sv_id]["End"][i])
            start2 = int(my_dict[sv_id]["Start_b"][i])
            end2 = int(my_dict[sv_id]["End_b"][i])
            sv_id2 = my_dict[sv_id]["SV_ID_b"][i]
            overlap = int(my_dict[sv_id]["Overlap"][i])
            fract1 = float(my_dict[sv_id]["Fraction"][i])
            fract2 = float(my_dict[sv_id]["Fraction_b"][i])

            if fract1 >= minf and fract2 >= minf:
                pass_filter[sv_id]["Start"].append(start)
                pass_filter[sv_id]["End"].append(end)
                pass_filter[sv_id]["Start_b"].append(start2)
                pass_filter[sv_id]["End_b"].append(end2)
                pass_filter[sv_id]["SV_ID_b"].append(sv_id2)
                pass_filter[sv_id]["Overlap"].append(overlap)
                pass_filter[sv_id]["Fraction"].append(fract1)
                pass_filter[sv_id]["Fraction_b"].append(fract2)

        if len(pass_filter[sv_id]["Start"]) == 0:
            del pass_filter[sv_id]

    return(pass_filter)

                
def get_best(my_dict,source,datas):
    ## determines "best" overlap via product of fraction and fraction_b
    ## whichever overlap has highest fraction product is returned
    ## in the event that multiple IDs are found to be "best"
    ## tries to run tiebreaker based on highest AF
    ## if tie persists, highest HomAlt count
    ## if tie still persists, returns the first ID

    best = defaultdict(list)
    final = defaultdict(list)

    for sv_id in my_dict:

        ids = []
        fracts = []
        for i in range(len(my_dict[sv_id]["Overlap"])):
            sv_id2 = my_dict[sv_id]["SV_ID_b"][i]
            ids.append(sv_id2)
            fract1 = my_dict[sv_id]["Fraction"][i]
            fract2 = my_dict[sv_id]["Fraction_b"][i]
            prod_fract = fract1 * fract2
            fracts.append(prod_fract)
            
        maxf = max(fracts)
        indices = [i for i, x in enumerate(fracts) if x == maxf]

        if len(indices) == 1:
            best[sv_id].append(str(ids[indices[0]]))

        elif len(indices) > 1:
            afs = []
            homalts = []
            tiebreaker_ids = []

            for i in indices:
                af = datas[source][str(ids[i])][0]
                afs.append(af)
                ha = datas[source][str(ids[i])][3]
                homalts.append(ha)
                tiebreaker_ids.append(str(ids[i]))

            maxaf = max(afs)
            af_tiebreaker = [i for i, x in enumerate(afs) if x == maxaf]

            if len(af_tiebreaker) == 1:
                best[sv_id].append(str(tiebreaker_ids[af_tiebreaker[0]]))

            elif len(af_tiebreaker) > 1:
                maxha = max(homalts)
                ha_tiebreaker = [i for i, x in enumerate(homalts) if x == maxha]

                if len(ha_tiebreaker) == 1:
                    best[sv_id].append(str(tiebreaker_ids[ha_tiebreaker[0]]))

                if len(ha_tiebreaker) > 1:
                    best[sv_id].append(str(ids[indices[0]]))
    return(best)


def create_max_annotations(my_list,vcf,req_sources):
    ## for each item in the list
    ## create Max_AF,Max_Het,Max_HomAlt annotations
    for i in my_list:
        af = ["Max",str(i),"AF"]
        vcf.add_info_to_header({"ID": "_".join(af), "Description": "The maximum " + str(i) + " AF from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Float", "Number": "1"})
        het = ["Max",str(i),"Het"]
        vcf.add_info_to_header({"ID": "_".join(het), "Description": "The maximum " + str(i) + " Het count from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Integer", "Number": "1"})
        homalt = ["Max",str(i),"HomAlt"]
        vcf.add_info_to_header({"ID": "_".join(homalt), "Description": "The maximum " + str(i) + " HomAlt count from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Integer", "Number": "1"})


def create_best_annotations(source,my_list,vcf):
    ## for each item in the list
    ## create Best_source_AF,Best_source_Het,Best_source_HomAlt annotations
    for i in my_list:
        af = ["Best",str(source),str(i),"AF"]
        vcf.add_info_to_header({"ID": "_".join(af), "Description": "The best " + str(i) + " AF match for " + str(source), "Type": "Float", "Number": "1"})
        het = ["Best",str(source),str(i),"Het"]
        vcf.add_info_to_header({"ID": "_".join(het), "Description": "The best " + str(i) + " Het count match for " + str(source), "Type": "Integer", "Number": "1"})
        homalt = ["Best",str(source),str(i),"HomAlt"]
        vcf.add_info_to_header({"ID": "_".join(homalt), "Description": "The best " + str(i) + " HomAlt count match " + str(source), "Type": "Integer", "Number": "1"})


#def join_annotations(source):
    ## create annotations specific to sources requested
    ## makes "best" annotations corresponding to specified fields
#    vcf.add_info_to_header({"ID": "Best_" + source, "Description": "The " + source + " ID of the best matching SV", "Type": "String", "Number": "1"})
#    for i in source_cols[source]:
#        annotation = bed_headers[i]
#        annotation_type = header_types[annotation]
#        vcf.add_info_to_header({"ID": "Best_" + source + "_" + annotation, "Description": "The best match " + annotation + " for " + source, "Type": annotation_type, "Number": "1"})


def get_feature(source,my_list,col,datas):
    ## takes a list of SV_ID_b and returns the feature for the specified col
    ## returns list of requested feature
    features = []
    for i in my_list:
        feature = datas[source][i][col]
        if feature != "NA": 
            features.append(feature)
    return(features)


#def add_best_annotations(source,my_list):
    ## writes out best annotations to vcf
#    if len(my_list) > 1:
#        print(source + " has too many best")
#    else:
#        v.INFO["Best_" + source] = my_list[0]
#        for i in source_cols[source]:
#            annotation = bed_headers[i]
#            annotation_type = header_types[annotation]
#            v.INFO["Best_" + source + "_" + annotation] = datas[source][my_list[0]][i]


def write_max_values(sv_id,my_dict,features,req_sources,datas,header_cols,header_types,v):
    ## features is a list of annotations
    ## finds max values for list of annotations
    ## adds those max values to output VCF
    for i in features:
        anno = "Max_" + str(i)
        max_val = float(0)
        all_vals = []
        for source in req_sources:
            matches_vals = []
            if sv_id in my_dict[source]:
                matches_vals = get_feature(source,my_dict[source][sv_id],header_cols[i],datas)
            all_vals.extend(matches_vals)
        if len(all_vals) > 0:
            if header_types[i] == "Float":
                max_val = max(list(map(float,all_vals)))
            if header_types[i] == "Integer":
                max_val = max(list(map(int,all_vals)))
        v.INFO[anno] = max_val


def write_best_values(source,sv_id,my_dict,features,datas,header_cols,v):
    ## features is a list of annotations
    ## finds best values for list of annotations
    ## adds those best values to output VCF
    for i in features:
        anno = "Best_" + str(source) + "_" + str(i)
        best_val = []
        if sv_id in my_dict[source]:
            matches_vals = get_feature(source,my_dict[source][sv_id],header_cols[i],datas)
            best_val.extend(matches_vals)
        if len(best_val) == 0:
            best_val = [0]
        v.INFO[anno] = best_val[0]


def create_full_string(source,my_list,datas):
    ## for a list of SV_ID_b
    ## generate comma-separated list of all annotations with |"s
    big_list = []
    for i in my_list:
        anno_list = [i]
        anno_list.extend(datas[source][i])
        big_string = "|".join(list(map(str,anno_list)))
        big_list.append(big_string)
    return(",".join(big_list))


def annotate(parser,args):

    ########################################
    ##        Interpret Arguments         ##
    ########################################

    ## Check number of CPUs
    ncpus = int(args.cpu) if int(args.cpu) > 0 else 1
    if ncpus > int(multiprocessing.cpu_count()): 
        print("\n**Warning** Too many CPUs designated. Number of CPUs will be set to {}".format( multiprocessing.cpu_count() if multiprocessing.cpu_count() > 0 else 1)) 
        ncpus = int(multiprocessing.cpu_count()) if int(multiprocessing.cpu_count())  > 0 else 1

    ## Check for source data file
    if args.pickled_source is None and args.bed is None:
        raise argparse.ArgumentTypeError("The following arguments are required: -b/--bed OR --pickled-source")    

    ## identify and set up input and output vcfs
    vcf = VCF(args.vcf, threads=ncpus)
    print("\nAnnotating the following VCF: {}".format(args.vcf))

    ## output vcf file
    output_vcf = args.out

    ## save specified minimum overlap fraction threshold
    ## if none provided use 0.001
    if args.minf is not None:
        minf = args.minf
        for f in minf:
            if float(f) > 1 or float(f) < 0:
                raise NameError("\nThe value included with -f, " + f + ", must be between 0.0 and 1.0")
    else:
        minf = [float(0.001)]
        print("\nNo reciprocal overlap fraction indicated; using default 0.001")

    ## local vars
    sources = []
    datas = {}
    bed_lists = {}
    bed_headers = []

    ## Get info from bed file
    if args.bed:
        sources, datas, bed_lists, bed_headers = process_bed_source(args.bed)
    elif args.pickled_source:
        sources, datas, bed_lists, bed_headers = process_pickled_source(args.pickled_source)
    else:
        raise NameError("Data Source is missing")

    header_types, header_cols = get_header_types(bed_headers)

    ## prepare specified sources
    ## use all sources in bed if none are provided
    ## align sources with specified minf values
    print("\nAnnotating with the following sources and reciprocal overlap fractions:")
    req_sources = args.sources if args.sources is not None else sources
    minfs = defaultdict(float)

    ## Check for bad sources
    for s in req_sources:
        if s not in sources:
            raise NameError(str(s) + " is an unexpected source; acceptable sources include (case-sensitive): " + ",".join(sorted(sources)))
        minfs[s] = float(0.001)

    ## Prepare sources with minf
    if len(minf) == 1:
        for s in req_sources:
            minfs[s] = float(minf[0])
            print(s + "\t" + str(minf[0]))
    elif len(minf) > 1:
        if len(minf) != len(req_sources):
            raise NameError("Please list either 1 overlap fraction for all sources or the same number of overlap fractions (-f) as there are requested sources (-s)")
        for i in range(0,len(req_sources)):
            minfs[req_sources[i]] = float(minf[i])
            print(req_sources[i] + "\t" + minf[i])

    extras = args.ann if args.ann is not None else []

    ########################################
    ##          Process Input VCF         ##
    ########################################

    ## make lists of SV coordinates to be used with pyranges
    ## set up dicts for storing results
    print("\nGathering SV coordinates from VCF file: {}".format(args.vcf))
    lists = defaultdict(list)
    for v in vcf:

        ## VCF Info
        cipos = v.INFO.get("CIPOS95") if args.ci is not None else 0
        ciend = v.INFO.get("CIEND95") if args.ci is not None else 0 
        svtype = v.INFO.get("SVTYPE") if v.INFO.get("SVTYPE") is not None else None
        sv_id = v.ID

        ## Get variant info
        ## Check for BNDs 
        ## Check for confidence intervals for adjustments
        chrom = v.CHROM[3:] if v.CHROM.startswith("chr") else v.CHROM
        start = int(v.start) 
        start += int(cipos[0]) if args.ci == "out" else int(cipos[1]) if args.ci == "in" else 0
        end = int(v.INFO.get("END")) if svtype != "BND" else int(v.POS)
        end += int(ciend[1]) if args.ci == "out" else int(ciend[0]) if args.ci == "out" else 0
        if args.ci == "in" and end < start:
            end = start + 1

        ## Check for problematic coordinates
        if start > end:
            print("Problematic coordinates (start > end):")
            print(str(chrom) + "\t" + str(start) + "\t" + str(end))
            fix_start = start
            fix_end = end
            start = fix_end
            end = fix_start
            print("Flipped to:")
            print(str(chrom) + "\t"+ str(start) + "\t" + str(end))

        features = ["chrom","start","end","svtype","sv_id"]
        variables = [chrom,start,end,svtype,sv_id]
        for i in range(len(features)):
            lists[features[i]].append(variables[i])

    ## Close input VCF
    ## Reopen input for writing to output
    vcf.close(); vcf = VCF(args.vcf, threads=ncpus)

    ########################################
    ##       New Annotations To Add       ##
    ########################################

    ## these are the columns in the bed with the annotations
    ## some sources have limited annotations
    ## specify the ones you want to add here
    ## for now this needs to be hardcoded
    source_cols = {}
    source_cols["CCDG"] = [0]
    source_cols["CEPH"] = list(list(range(0,14,1)) + [98,99])
    source_cols["gnomAD"] = list(list(range(0,84,1)) + [98,99])
    source_cols["1000G_Smoove"] = list(list(range(0,70,1)) + list(range(84,100,1)))

    ## main, default max annotations to add
    vcf.add_info_to_header({"ID": "Max_AF", "Description": "The maximum AF from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Float", "Number": "1"})
    vcf.add_info_to_header({"ID": "Max_Het", "Description": "The maximum Het count from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Integer", "Number": "1"})
    vcf.add_info_to_header({"ID": "Max_HomAlt", "Description": "The maximum HomAlt count from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Integer", "Number": "1"})
    vcf.add_info_to_header({"ID": "Max_PopMax_AF", "Description": "The maximum PopMax_AF from all matching SVs across all specified data sources (" + ",".join(req_sources) + ")", "Type": "Float", "Number": "1"})

    ## add male and female annotations
    if "mf" in extras or "all" in extras:
        create_max_annotations(["Male","Female"],vcf,req_sources)

    ## add AFR annotations
    if "AFR" in extras or "pops" in extras or "all" in extras:
        create_max_annotations(["AFR"],vcf,req_sources)
        if "mf" in extras or "all" in extras:
            create_max_annotations(["AFR_Male","AFR_Female"],vcf,req_sources)

    ## add AMR annotations
    if "AMR" in extras or "pops" in extras or "all" in extras:
        create_max_annotations(["AMR"],vcf,req_sources)
        if "mf" in extras or "all" in extras:
            create_max_annotations(["AMR_Male","AMR_Female"],vcf,req_sources)

    ## add EAS annotations
    if "EAS" in extras or "pops" in extras or "all" in extras:
        create_max_annotations(["EAS"],vcf,req_sources)
        if "mf" in extras or "all" in extras:
            create_max_annotations(["EAS_Male","EAS_Female"],vcf,req_sources)

    ## add EUR annotations
    if "EUR" in extras or "pops" in extras or "all" in extras:
        create_max_annotations(["EUR"],vcf,req_sources)
        if "mf" in extras or "all" in extras:
            create_max_annotations(["EUR_Male","EUR_Female"],vcf,req_sources)

    ## add OTH annotations
    if "OTH" in extras or "pops" in extras or "all" in extras:
        create_max_annotations(["OTH"],vcf,req_sources)
        if "mf" in extras or "all" in extras:
            create_max_annotations(["OTH_Male","OTH_Female"],vcf,req_sources)

    ## add SAS annotations
    if "SAS" in extras or "pops" in extras or "all" in extras:
        create_max_annotations(["SAS"],vcf,req_sources)
        if "mf" in extras or "all" in extras:
            create_max_annotations(["SAS_Male","SAS_Female"],vcf,req_sources)

    for source in req_sources:
        vcf.add_info_to_header({"ID": source + "_Count", "Description": "The number of matching SVs with " + source, "Type": "Integer", "Number": "1"})

        ## add best annotations
        if "best" in extras or "all" in extras:
            vcf.add_info_to_header({"ID": "Best_" + source + "_ID", "Description": "The " + source + " ID of the best matching SV", "Type": "String", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_AF", "Description": "The best AF match for " + str(source), "Type": "Float", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_Het", "Description": "The best Het count match for " + str(source), "Type": "Integer", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_HomAlt", "Description": "The best HomAlt count match for " + str(source), "Type": "Integer", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_PopMax_AF", "Description": "The best PopMax_AF match for " + str(source), "Type": "Float", "Number": "1"})

            ## add male and female best annotations
            if "mf" in extras or "all" in extras:
                create_best_annotations(source,["Male","Female"],vcf)

            ## add AFR best annotations
            if "AFR" in extras or "pops" in extras or "all" in extras:
                create_best_annotations(source,["AFR"],vcf)
                if "mf" in extras or "all" in extras:
                    create_best_annotations(source,["AFR_Male","AFR_Female"],vcf)

            ## add AMR best annotations
            if "AMR" in extras or "pops" in extras or "all" in extras:
                create_best_annotations(source,["AMR"],vcf)
                if "mf" in extras or "all" in extras:
                    create_best_annotations(source,["AMR_Male","AMR_Female"],vcf)

            ## add EAS best annotations
            if "EAS" in extras or "pops" in extras or "all" in extras:
                create_best_annotations(source,["EAS"],vcf)
                if "mf" in extras or "all" in extras:
                    create_best_annotations(source,["EAS_Male","EAS_Female"],vcf)

            ## add EUR best annotations
            if "EUR" in extras or "pops" in extras or "all" in extras:
                create_best_annotations(source,["EUR"],vcf)
                if "mf" in extras or "all" in extras:
                    create_best_annotations(source,["EUR_Male","EUR_Female"],vcf)

            ## add OTH best annotations
            if "OTH" in extras or "pops" in extras or "all" in extras:
                create_best_annotations(source,["OTH"],vcf)
                if "mf" in extras or "all" in extras:
                    create_best_annotations(source,["OTH_Male","OTH_Female"],vcf)

            ## add SAS best annotations
            if "SAS" in extras or "pops" in extras or "all" in extras:
                create_best_annotations(source,["SAS"],vcf)
                if "mf" in extras or "all" in extras:
                    create_best_annotations(source,["SAS_Male","SAS_Female"],vcf)
        
        ## add mismatch annotations
        if "mis" in extras or "all" in extras:
            vcf.add_info_to_header({"ID": source + "_Mismatches", "Description": "Comma-separated list of the " + source + " IDs of overlapping SVs with different SVTYPEs", "Type": "String", "Number": "."})
            vcf.add_info_to_header({"ID": source + "_Mismatches_Count", "Description": "The number of  " + source + " overlapping SVs with different SVTYPEs", "Type": "Integer", "Number": "1"})
            vcf.add_info_to_header({"ID": source + "_Mismatch_SVTYPEs", "Description": "Comma-separated list of the other overlapping SVTYPEs for ", "Type": "String", "Number": "."})
            vcf.add_info_to_header({"ID": "Best_" + source + "_Mismatch_ID", "Description": "The " + source + " ID of the best overlapping SV with different SVTYPE", "Type": "String", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_Mismatch_SVTYPE", "Description": "The " + source + " SVTYPE of the best overlapping SV with different SVTYPE", "Type": "String", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_Mismatch_AF", "Description": "The " + source + " AF for the best overlapping SV with different SVTYPE", "Type": "Float", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_Mismatch_Het", "Description": "The " + source + " Het count for the best overlapping SV with different SVTYPE", "Type": "Integer", "Number": "1"})
            vcf.add_info_to_header({"ID": "Best_" + source + "_Mismatch_HomAlt", "Description": "The " + source + " HomAlt count for the best overlapping SV with different SVTYPE", "Type": "Integer", "Number": "1"})

        ## add full list of matches annotations
        if "full" in extras or "all" in extras:
            vcf.add_info_to_header({"ID": source + "_Matches", "Description": "Comma-separated list of  each SV match with " + source + ", where all annotations are listed as SV_ID|AF|HomRef|Het|HomAlt|Male_AF|Male_HomRef|Male_Het|Male_HomAlt|HemiAlt|Hemi_AF|Female_AF|Female_HomRef|Female_Het|Female_HomAlt|AFR_AF|AFR_HomRef|AFR_Het|AFR_HomAlt|AFR_Male_AF|AFR_Male_HomRef|AFR_Male_Het|AFR_Male_HomAlt|AFR_HemiAlt|AFR_Hemi_AF|AFR_Female_AF|AFR_Female_HomRef|AFR_Female_Het|AFR_Female_HomAlt|AMR_AF|AMR_HomRef|AMR_Het|AMR_HomAlt|AMR_Male_AF|AMR_Male_HomRef|AMR_Male_Het|AMR_Male_HomAlt|AMR_HemiAlt|AMR_Hemi_AF|AMR_Female_AF|AMR_Female_HomRef|AMR_Female_Het|AMR_Female_HomAlt|EAS_AF|EAS_HomRef|EAS_Het|EAS_HomAlt|EAS_Male_AF|EAS_Male_HomRef|EAS_Male_Het|EAS_Male_HomAlt|EAS_HemiAlt|EAS_Hemi_AF|EAS_Female_AF|EAS_Female_HomRef|EAS_Female_Het|EAS_Female_HomAlt|EUR_AF|EUR_HomRef|EUR_Het|EUR_HomAlt|EUR_Male_AF|EUR_Male_HomRef|EUR_Male_Het|EUR_Male_HomAlt|EUR_HemiAlt|EUR_Hemi_AF|EUR_Female_AF|EUR_Female_HomRef|EUR_Female_Het|EUR_Female_HomAlt|OTH_AF|OTH_HomRef|OTH_Het|OTH_HomAlt|OTH_Male_AF|OTH_Male_HomRef|OTH_Male_Het|OTH_Male_HomAlt|OTH_HemiAlt|OTH_Hemi_AF|OTH_Female_AF|OTH_Female_HomRef|OTH_Female_Het|OTH_Female_HomAlt|SAS_AF|SAS_HomRef|SAS_Het|SAS_HomAlt|SAS_Male_AF|SAS_Male_HomRef|SAS_Male_Het|SAS_Male_HomAlt|SAS_HemiAlt|SAS_Hemi_AF|SAS_Female_AF|SAS_Female_HomRef|SAS_Female_Het|SAS_Female_HomAlt|PopMax_AF|In_Pop as found in " + args.bed if args.bed else args.pickled_source, "Type": "String", "Number": "."})


    #if "CCDG" in req_sources:
    #    join_annotations("CCDG")
    #if "CEPH" in req_sources:
    #    join_annotations("CEPH")
    #if "gnomAD" in req_sources:
    #    join_annotations("gnomAD")
    #if "1000G_Smoove" in req_sources:
    #    join_annotations("1000G_Smoove")

    ########################################
    ##        Run pyRanges Commands       ##
    ########################################

    print("\nCreating pyranges object from VCF data...")
    pr_vcf = pr.from_dict({"Chromosome": lists["chrom"], \
                           "Start": lists["start"], \
                           "End": lists["end"], \
                           "SVTYPE": lists["svtype"], \
                           "SV_ID": lists["sv_id"]})

    ## pyRanges join
    ## identify intersections of SV coordinates
    join_matches = {}
    join_best_matches = {}
    join_mismatches = {}
    join_best_mismatches = {}

    for source in req_sources:
        print("\nRunning pyranges.join for:")
        print(source)
        pr_bed = pr.from_dict({"Chromosome": bed_lists[source]["chrom"], \
                           "Start": bed_lists[source]["start"], \
                           "End": bed_lists[source]["end"], \
                           "SVTYPE": bed_lists[source]["svtype"], \
                           "SV_ID": bed_lists[source]["sv_id"]})
        joined = pr_vcf.join(pr_bed, how = "left").sort()

        pr_matches = joined.subset(lambda df: (df.SVTYPE == df.SVTYPE_b))
        pr_mismatches = joined.subset(lambda df: (df.SVTYPE != df.SVTYPE_b)).subset(lambda df: df.SVTYPE_b != "-1")
        pr_noBND = joined.subset(lambda df: (df.SVTYPE != df.SVTYPE_b)).subset(lambda df: df.SVTYPE_b != "-1").subset(lambda df: df.SVTYPE_b != "BND")
        pr_uniqs = joined.subset(lambda df: df.SVTYPE_b == "-1")

        ##NOTE: Using mutliple CPUs for pyranges slows down the process 

        ## create SVTYPE matches dict
        if source not in join_matches:
            join_matches[source] = {}
        matches = convert_dict(pr_matches)
        filtered_matches = reciprocal_overlap(matches,source,minfs)
        filtered_matches_ids = {}
        for sv_id in filtered_matches:
            if sv_id not in filtered_matches_ids:
                filtered_matches_ids[sv_id] = []
            for i in filtered_matches[sv_id]["SV_ID_b"]:
                filtered_matches_ids[sv_id].append(i)
            if sv_id not in join_matches[source]:
                join_matches[source][sv_id] = []
            join_matches[source][sv_id].extend(filtered_matches_ids[sv_id])

        ## create SVTYPE best match dict
        if source not in join_best_matches:
            join_best_matches[source] = {}
        best_matches_ids = get_best(filtered_matches,source,datas)
        for sv_id in best_matches_ids:
            if sv_id not in join_best_matches[source]:
                join_best_matches[source][sv_id] = []
            join_best_matches[source][sv_id].extend(best_matches_ids[sv_id])

        ## create SVTYPE mismatches dict
        if source not in join_mismatches:
            join_mismatches[source] = {}
        mismatches = convert_dict(pr_mismatches)
        filtered_mismatches = reciprocal_overlap(mismatches,source,minfs)
        filtered_mismatches_ids = {}
        for sv_id in filtered_mismatches:
            if sv_id not in filtered_mismatches_ids:
                filtered_mismatches_ids[sv_id] = []
            for i in filtered_mismatches[sv_id]["SV_ID_b"]:
                filtered_mismatches_ids[sv_id].append(i)
            if sv_id not in join_mismatches[source]:
                join_mismatches[source][sv_id] = []
            join_mismatches[source][sv_id].extend(filtered_mismatches_ids[sv_id])

        ## create SVTYPE best mismatch dict
        if source not in join_best_mismatches:
            join_best_mismatches[source] = {}
        best_mismatches_ids = get_best(filtered_mismatches,source,datas)
        for sv_id in best_mismatches_ids:
            if sv_id not in join_best_mismatches[source]:
                join_best_mismatches[source][sv_id] = []
            join_best_mismatches[source][sv_id].extend(best_mismatches_ids[sv_id])

    ########################################
    ## Write To Output/Adding Annotations ##
    ########################################

    print("\nAdding annotations to output VCF...")
    new_vcf = Writer(output_vcf, vcf)
    for v in vcf:
        sv_id = v.ID
        write_max_values(sv_id,join_matches,["AF","Het","HomAlt","PopMax_AF"],req_sources,datas,header_cols,header_types,v)

        ## write male female annotations 
        if "mf" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["Male_AF","Male_Het","Male_HomAlt","Female_AF","Female_Het","Female_HomAlt"],req_sources,datas,header_cols,header_types,v)

        ## write AFR annotations
        if "AFR" in extras or "pops" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["AFR_AF","AFR_Het","AFR_HomAlt"],req_sources,datas,header_cols,header_types,v)
            if "mf" in extras or "all" in extras:
                write_max_values(sv_id,join_matches,["AFR_Male_AF","AFR_Male_Het","AFR_Male_HomAlt","AFR_Female_AF","AFR_Female_Het","AFR_Female_HomAlt"],req_sources,datas,header_cols,header_types,v)

        ## write AMR annotations
        if "AMR" in extras or "pops" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["AMR_AF","AMR_Het","AMR_HomAlt"],req_sources,datas,header_cols,header_types,v)
            if "mf" in extras or "all" in extras:
                write_max_values(sv_id,join_matches,["AMR_Male_AF","AMR_Male_Het","AMR_Male_HomAlt","AMR_Female_AF","AMR_Female_Het","AMR_Female_HomAlt"],req_sources,datas,header_cols,header_types,v)

        ## write EAS annotations
        if "EAS" in extras or "pops" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["EAS_AF","EAS_Het","EAS_HomAlt"],req_sources,datas,header_cols,header_types,v)
            if "mf" in extras or "all" in extras:
                write_max_values(sv_id,join_matches,["EAS_Male_AF","EAS_Male_Het","EAS_Male_HomAlt","EAS_Female_AF","EAS_Female_Het","EAS_Female_HomAlt"],req_sources,datas,header_cols,header_types,v)

        ## write EUR annotations
        if "EUR" in extras or "pops" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["EUR_AF","EUR_Het","EUR_HomAlt"],req_sources,datas,header_cols,header_types,v)
            if "mf" in extras or "all" in extras:
                write_max_values(sv_id,join_matches,["EUR_Male_AF","EUR_Male_Het","EUR_Male_HomAlt","EUR_Female_AF","EUR_Female_Het","EUR_Female_HomAlt"],req_sources,datas,header_cols,header_types,v)

        ## write OTH annotations
        if "OTH" in extras or "pops" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["OTH_AF","OTH_Het","OTH_HomAlt"],req_sources,datas,header_cols,header_types,v)
            if "mf" in extras or "all" in extras:
                write_max_values(sv_id,join_matches,["OTH_Male_AF","OTH_Male_Het","OTH_Male_HomAlt","OTH_Female_AF","OTH_Female_Het","OTH_Female_HomAlt"],req_sources,datas,header_cols,header_types,v)

        ## write SAS annotations
        if "SAS" in extras or "pops" in extras or "all" in extras:
            write_max_values(sv_id,join_matches,["SAS_AF","SAS_Het","SAS_HomAlt"],req_sources,datas,header_cols,header_types,v)
            if "mf" in extras or "all" in extras:
                write_max_values(sv_id,join_matches,["SAS_Male_AF","SAS_Male_Het","SAS_Male_HomAlt","SAS_Female_AF","SAS_Female_Het","SAS_Female_HomAlt"],req_sources,datas,header_cols,header_types,v)
        
        counts = {}
        for source in req_sources:
        
            ## write data source matches count
            if source not in counts:
                counts[source] = 0
            if sv_id in join_matches[source]:
                counts[source] = len(join_matches[source][sv_id])
            v.INFO[source + "_Count"] = counts[source]
            
            ## write best annotations
            if "best" in extras or "all" in extras:
                if sv_id in join_best_matches[source]:
                    annoID = "Best_" + source + "_ID"
                    v.INFO[annoID] = join_best_matches[source][sv_id][0]
                    write_best_values(source,sv_id,join_best_matches,["AF","Het","HomAlt","PopMax_AF"],datas,header_cols,v)

                    ## write best male female annotations 
                    if "mf" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["Male_AF","Male_Het","Male_HomAlt","Female_AF","Female_Het","Female_HomAlt"],datas,header_cols,v)

                    ## write AFR best annotations
                    if "AFR" in extras or "pops" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["AFR_AF","AFR_Het","AFR_HomAlt"],datas,header_cols,v)
                        if "mf" in extras or "all" in extras:
                            write_best_values(source,sv_id,join_best_matches,["AFR_Male_AF","AFR_Male_Het","AFR_Male_HomAlt","AFR_Female_AF","AFR_Female_Het","AFR_Female_HomAlt"],datas,header_cols,v)

                    ## write AMR best annotations
                    if "AMR" in extras or "pops" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["AMR_AF","AMR_Het","AMR_HomAlt"],datas,header_cols,v)
                        if "mf" in extras or "all" in extras:
                            write_best_values(source,sv_id,join_best_matches,["AMR_Male_AF","AMR_Male_Het","AMR_Male_HomAlt","AMR_Female_AF","AMR_Female_Het","AMR_Female_HomAlt"],datas,header_cols,v)

                    ## write EAS best annotations
                    if "EAS" in extras or "pops" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["EAS_AF","EAS_Het","EAS_HomAlt"],datas,header_cols,v)
                        if "mf" in extras or "all" in extras:
                            write_best_values(source,sv_id,join_best_matches,["EAS_Male_AF","EAS_Male_Het","EAS_Male_HomAlt","EAS_Female_AF","EAS_Female_Het","EAS_Female_HomAlt"],datas,header_cols,v)

                    ## write EUR best annotations
                    if "EUR" in extras or "pops" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["EUR_AF","EUR_Het","EUR_HomAlt"],datas,header_cols,v)
                        if "mf" in extras or "all" in extras:
                            write_best_values(source,sv_id,join_best_matches,["EUR_Male_AF","EUR_Male_Het","EUR_Male_HomAlt","EUR_Female_AF","EUR_Female_Het","EUR_Female_HomAlt"],datas,header_cols,v)

                    ## write OTH best annotations
                    if "OTH" in extras or "pops" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["OTH_AF","OTH_Het","OTH_HomAlt"],datas,header_cols,v)
                        if "mf" in extras or "all" in extras:
                            write_best_values(source,sv_id,join_best_matches,["OTH_Male_AF","OTH_Male_Het","OTH_Male_HomAlt","OTH_Female_AF","OTH_Female_Het","OTH_Female_HomAlt"],datas,header_cols,v)

                    ## write SAS best annotations
                    if "SAS" in extras or "pops" in extras or "all" in extras:
                        write_best_values(source,sv_id,join_best_matches,["SAS_AF","SAS_Het","SAS_HomAlt"],datas,header_cols,v)
                        if "mf" in extras or "all" in extras:
                            write_best_values(source,sv_id,join_best_matches,["SAS_Male_AF","SAS_Male_Het","SAS_Male_HomAlt","SAS_Female_AF","SAS_Female_Het","SAS_Female_HomAlt"],datas,header_cols,v)

            ## write mismatch annotations
            if "mis" in extras or "all" in extras:
                if sv_id in join_mismatches[source]:
                    v.INFO[source + "_Mismatches"] = ",".join(join_mismatches[source][sv_id])
                    v.INFO[source + "_Mismatches_Count"] = len(join_mismatches[source][sv_id])
                    mismatch_svtypes = get_feature(source,join_mismatches[source][sv_id],header_cols["SVTYPE"],datas)
                    v.INFO[source + "_Mismatch_SVTYPEs"] = ",".join(sorted(list(set(mismatch_svtypes))))
                    v.INFO["Best_" + source + "_Mismatch_ID"] = join_best_mismatches[source][sv_id][0]
                    mismatch_svtype = get_feature(source,join_best_mismatches[source][sv_id],header_cols["SVTYPE"],datas)
                    v.INFO["Best_" + source + "_Mismatch_SVTYPE"] = mismatch_svtype[0]
                    mismatch_AF = get_feature(source,join_best_mismatches[source][sv_id],header_cols["AF"],datas)
                    v.INFO["Best_" + source + "_Mismatch_AF"] = mismatch_AF[0]
                    mismatch_Het = get_feature(source,join_best_mismatches[source][sv_id],header_cols["Het"],datas)
                    if len(mismatch_Het) == 0:
                        mismatch_Het = [0]
                    v.INFO["Best_" + source + "_Mismatch_Het"] = mismatch_Het[0]
                    mismatch_HomAlt = get_feature(source,join_best_mismatches[source][sv_id],header_cols["HomAlt"],datas)
                    if len(mismatch_HomAlt) == 0:
                        mismatch_HomAlt = [0]
                    v.INFO["Best_" + source + "_Mismatch_HomAlt"] = mismatch_HomAlt[0]

            ## write full matches annotations
            if "full" in extras or "all" in extras:
                if sv_id in join_matches[source]:
                    v.INFO[source + "_Matches"] = create_full_string(source,join_matches[source][sv_id],datas)            
                    
        new_vcf.write_record(v)

    new_vcf.close(); vcf.close()

    print("\nAnnoted VCF writen to: {}".format(args.out))
    print("\nDONE")

