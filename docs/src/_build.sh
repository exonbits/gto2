#!/bin/sh

set -ev

Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"

rm -r ../book/libs/ ../book/* || echo "nothing to remove"
rm -r _book/book/images || echo "nothing to remove"
mv _book/* ../book/
rmdir _book
