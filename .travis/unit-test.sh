#!/bin/bash

set -e
errors=0

# Run unit tests
python MC_Star/MC_Star_test.py || {
    echo "'python python/MC_Star/MC_Star_test.py' failed"
    let errors+=1
}

# Check program style
pylint -E MC_Star/*.py || {
    echo 'pylint -E MC_Star/*.py failed'
    let errors+=1
}

[ "$errors" -gt 0 ] && {
    echo "There were $errors errors found"
    exit 1
}

echo "Ok : Python specific tests"
