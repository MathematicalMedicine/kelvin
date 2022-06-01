#!/bin/bash

githash=$(git describe --always 2>/dev/null)$(git status 2>/dev/null | grep "modified:\|added:\|deleted:" -q && echo "-M")
if [[ ${githash-} ]]; then
    echo $githash
else
    echo "(none)"
fi
