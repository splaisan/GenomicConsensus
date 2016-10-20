#!/bin/bash
set -euo pipefail

echo "# DEPENDENCIES"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load git gcc/5.3.0 python/2.7.9 cmake cram swig ccache virtualenv zlib/1.2.5 ninja boost

set +u
source venv/bin/activate
set -u

echo "# TEST"

echo "## CC2 version test"
python -c "import ConsensusCore2 ; print ConsensusCore2.__version__"

echo "## test GC"
make check

echo "## test GC extra"
make extra-tests

echo "## test GC internal"
make internal-tests

set +u
deactivate
set -u
