#!/bin/sh
MATMUL_DIR=`pwd`

cat > serial.qsub <<EOF
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

export LD_LIBRARY_PATH=/share/apps/local/lib64:\$LD_LIBRARY_PATH
cd $MATMUL_DIR
hostname > timing-$2.out
$1 >> timing-$2.out
exit 0;
EOF
