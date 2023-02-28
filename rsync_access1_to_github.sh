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
git add README.txt ;
git add dir_m/*.m ;
git add dir_m_dependencies/*.m ;
git add *.sh ;
git commit -m "updating dir_PAD with test_Chimpanzee_4 " ;
git push ;
git pull ;
cd /data/rangan/dir_bcc/dir_PAD ;

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
#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
#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
#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
#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
#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
#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
#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
