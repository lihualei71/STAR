#!/bin/bash

params_filename='corr_params_convex.txt'

while read rho type times seed
do
    filename='../results/STAR_convex_rho_"'$rho'"_type_"'$type'"_times_"'$times'"_seed_"'$seed'".sh'
    echo $filename
    Rout_filename='../results/STAR_convex_rho_"'$rho'"_type_"'$type'"__times_"'$times'"_seed_"'$seed'".Rout'
    echo $Rout_filename
    touch $filename
    echo '#!/bin/bash' > $filename
    echo '' >> $filename
    echo 'export rho="'$rho'"' >> $filename
    echo 'export type="'$type'"' >> $filename
    echo 'export times="'$times'"' >> $filename
    echo 'export seed="'$seed'"' >> $filename    
    echo "R CMD BATCH --no-save STAR_convex_expr.R "$Rout_filename >> $filename
    chmod 755 $filename
    sbatch -p high $filename
done < $params_filename
