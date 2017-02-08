echo "#pragma once" > glyphs.h
echo "#define FONT_WIDTH 7" >> glyphs.h
echo "#define FONT_HEIGHT 12" >> glyphs.h
echo "static unsigned char font_digits[][12] = {" >> glyphs.h
for cid in {33..127}
do
	if [ "$cid" == 127 ]
	then
		C=`echo "°" | iconv -t iso-8859-1`
		PARAM=
		CUT="-left 14 -top 10"
	else
		C=$(printf "\x$(printf '%x' $cid)")
		PARAM="-builtin fixed"
		CUT="-top 6 -left 8"
	fi
	STR="$STR$C";
	pbmtext $PARAM "$C" | pnmcut $CUT -height 12 -width 7 >../temp/temp.pbm 
	convert ../temp/temp.pbm ../temp/temp.xbm
	echo "{ $(grep 0x ../temp/temp.xbm) },"  >> glyphs.h
done
echo "};" >> glyphs.h

rm -f ../temp/temp.pbm ../temp/temp.xbm

