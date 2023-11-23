input_folder = os.environ['INPUT_DIR']
output_folder = os.environ['OUTPUT_DIR']
env_folder = os.environ['ENV_DIR']
script_folder = os.environ['SCRIPTS_DIR']
submodule_folder = os.environ['SUBMODULES_DIR']

map_labels = ["EarlyDryPoor", "EarlyWetRich", "EarlyWetPoor", "LateDryPoor", "EarlyDryRich", "LateWetPoor", "LateDryRich", "LateWetRich"]
net_types = ["jaccard"]

include: "rules/parser.rule"
include: "rules/jaccard.rule"

rule all:
    input:
        #output_folder + "/PA_combined_table/PA_OTUs.csv",
        output_folder + "/PA_split_nets/combined_nets.xlsx",


