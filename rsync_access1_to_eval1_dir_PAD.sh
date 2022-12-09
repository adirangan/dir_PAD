rsync -avum \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/dir_PAD \
    /home/rangan/dir_bcc/ \
    --include="*/" \
    --include="*.h" \
    --include="*.c" \
    --include="*.in" \
    --include="*.make" \
    --include="*.m" \
    --include="*.sh" \
    --exclude="*" ;
