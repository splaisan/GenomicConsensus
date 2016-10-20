#!/bin/bash

if [ "$(ps -p "$$" -o comm=)" != "bash" ]; then
    # Taken from http://unix-linux.questionfor.info/q_unix-linux-programming_85038.html
    bash "$0" "$@"
    exit "$?"
fi

echo "# DEPENDENCIES"
echo "## Load modules"
source /mnt/software/Modules/current/init/bash
module load git gcc/5.3.0 python/2.7.9 cmake cram swig ccache virtualenv zlib/1.2.5 ninja boost

source venv/bin/activate

echo "# TEST"

echo "## CC2 version test"
python -c "import ConsensusCore2 ; print ConsensusCore2.__version__"

echo "## test GC"
make check

echo "## test GC extra"
make extra-tests

echo "## test GC internal"
make internal-tests

deactivate
