rsync -avum \
      /data/rangan/dir_bcc/dir_PAD \
      /home/rangan/dir_bcc/ \
      --include="*/" \
      --include="*.sh" \
      --include="*.m" \
      --include="*.c" \
      --include="*.h" \
      --include="*.f" \
      --include="*.tex" \
      --include="*.bib" \
      --include="*.pdf" \
      --include="*.fig" \
      --exclude="*" \
      ;
