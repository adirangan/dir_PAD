rsync -avum \
    /home/rangan/dir_bcc/dir_PAD \
    rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/ \
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
# 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
# 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
