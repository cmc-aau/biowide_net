# biowide_net

Repo for the creation and analysis of the ecological network for the [Biowide](https://ecos.au.dk/en/researchconsultancy/themes/biowide) project.

Please copy the input data in the `/data` folder and amend the `/config/main.yaml` file accordingly (if required).

The script requires mamba and conda but installs all the other required packages via mamba. The results from the paper can be reproduced by running:
```
bash scripts/scripts_bash/PA_net_analysis.sh
```

The results will be collected in the `/results` folder.

