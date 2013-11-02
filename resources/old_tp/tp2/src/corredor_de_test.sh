make new
#echo -e "ejemplo1\n"
#./main casos_test/ejemplo1-pagelinks.tp2 casos_test/ejemplo1-page.tp2
#echo -e "\n"
echo -e "ejemplo2\n"
./main -directo casos_test/ejemplo2-pagelinks.tp2 casos_test/ejemplo2-page.tp2 ranking.tp2
#echo -e "\n"
#echo -e "dc4 \n"
#./main casos_test/dc4-pagelinks.tp2 casos_test/dc4-page.tp2
#echo -e "\n"
#echo -e "dc9 \n"
#./main casos_test/dc9-pagelinks.tp2 casos_test/dc9-page.tp2
#echo -e "\n"
#echo -e "dc22 \n"
#./main casos_test/dc22-pagelinks.tp2 casos_test/dc22-page.tp2

./main -iterativo casos_test/extwiki-20100920-pagelinks.tp2 casos_test/extwiki-20100920-page.tp2 ranking_gl.tp2 5000 0.00000001 >> resultado_iteraciones

#09-2010-text-www.unesco.org-pagelinks.tp2  dc4-page.tp2                   eswiki-20100912-page.tp2                        obama.transition.05.11.2009-text-page.tp2
#09-2010-text-www.unesco.org-page.tp2       dc9-pagelinks.tp2              extwiki-20100920-pagelinks.tp2                  simplewiki-20100902-pagelinks.tp2
#2008-Election-10-28-text-pagelinks.tp2     dc9-page.tp2                   extwiki-20100920-page.tp2                       simplewiki-20100902-page.tp2
#2008-Election-10-28-text-page.tp2          ejemplo1-pagelinks.tp2         glwiki-20100912-pagelinks.tp2                   US.CS.2004-11-text-pagelinks.tp2
#archivines.txt                             ejemplo1-page.tp2              glwiki-20100912-page.tp2                        US.CS.2004-11-text-page.tp2
#dc22-pagelinks.tp2                         ejemplo2-pagelinks.tp2         iawiki-20100910-pagelinks.tp2
#dc22-page.tp2                              ejemplo2-page.tp2              iawiki-20100910-page.tp2
#dc4-pagelinks.tp2                          eswiki-20100912-pagelinks.tp2  obama.transition.05.11.2009-text-pagelinks.tp2
