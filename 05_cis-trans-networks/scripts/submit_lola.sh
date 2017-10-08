#!/bin/bash

jid1=$(sbatch lola_cpg_core.sh)
sbatch lola_cpg_ext.sh
sbatch --dependency=afterok:$jid1 lola_snp_core.sh
sbatch --dependency=afterok:$jid1 lola_snp_ext.sh
