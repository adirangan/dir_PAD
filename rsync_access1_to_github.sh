rsync -avum \
    /data/rangan/dir_bcc/dir_PAD \
    /home/rangan/dir_bcc/ \
    --include="*/" \
    --include="*.h" \
    --include="*.c" \
    --include="*.in" \
    --include="*.make" \
    --include="*.m" \
    --include="*.sh" \
    --exclude="*" ;
#rsync -avum \
#    /home/rangan/dir_bcc/dir_PAD \
#    rangan@access1.cims.nyu.edu:/data/rangan/dir_bcc/ \
#    --include="*/" \
#    --include="*.h" \
#    --include="*.c" \
#    --include="*.in" \
#    --include="*.make" \
#    --include="*.m" \
#    --include="*.sh" \
#    --exclude="*" ;
cd /data/rangan/dir_bcc/dir_PAD ;
git add dir_m/*.m ;
git add dir_m_dependencies/*.m ;
git add *.sh ;
git commit -m "updating dir_PAD for shreya thirumalai " ;
git push ;
git pull ;
cd /home/rangan/dir_bcc/dir_PAD ;

#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
#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
#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
#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
#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
#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
#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
#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
#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
#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
#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
