#!/usr/bin/awk -f
BEGIN {
    RS = "(\r\n|\n\r|\r|\n)"
    FS = " *\t *"
    SUBSEP = ":"
}
FNR==1 {
    ++file
}
NF>=2 {
    if ($1 in keynum)
        key = keynum[$1]
    else {
        key = ++keys
        keynum[$1] = key
        keystr[key] = $1
    }
    printf "key = %s, file = %s, value = %s\n", key, file, $2 > "/dev/null"
    value[key,file] = $2
}
END {
    files = file
    for (key = 1; key <= keys; key++) {
        printf "%s", keystr[key]
        for (file = 1; file <= files; file++)
            printf "\t%s", value[key,file]
        printf "\n"
    }
}
