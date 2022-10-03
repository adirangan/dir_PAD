rsync -avum \
    /home/rangan/dir_bcc/dir_PAD \
    /data/rangan/dir_bcc/ \
    --include="*/" \
    --include="*.h" \
    --include="*.c" \
    --include="*.in" \
    --include="*.make" \
    --include="*.m" \
    --include="*.sh" \
    --exclude="*" ;
cd /home/rangan/dir_bcc/dir_PAD ;
git add dir_m/*.m ;
git add dir_m_dependencies/*.m ;
git add *.sh ;
git commit -m "updating dir_PAD for shreya thirumalai " ;
git push ;
git pull ;
cd /home/rangan/dir_bcc/dir_PAD ;

#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
