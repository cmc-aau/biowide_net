input_folder = "/home/fdelogu/biowide_net/data"
output_folder = "/home/fdelogu/biowide_net/results"
env_folder = "/home/fdelogu/biowide_net/envs"
script_folder = "/home/fdelogu/biowide_net/scripts"
submodule_folder = "/home/fdelogu/biowide_net/submodules"

map_labels = ["EarlyDryPoor", "EarlyWetRich", "EarlyWetPoor", "LateDryPoor", "EarlyDryRich", "LateWetPoor", "LateDryRich", "LateWetRich"]
net_types = ["sparcc", "jaccard", "pearson"]

include: "rules/parser.rule"
include: "rules/sparcc.rule"
include: "rules/jaccard.rule"
include: "rules/pearson.rule"

rule all:
    input:
        #expand(output_folder + "/split_nets/{net_type}/{map_label}_final.xlsx", map_label = map_labels, net_type = net_types)
        output_folder + "/split_nets/combined_nets.xlsx"


