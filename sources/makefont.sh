GLYPHS=glyphs.h
echo "#pragma once" > $GLYPHS
echo "#define FONT_WIDTH 16" >> $GLYPHS
echo "#define FONT_HEIGHT 16" >> $GLYPHS
echo "static unsigned char font_digits[][FONT_HEIGHT * (FONT_WIDTH + 7) / 8] = {" >> $GLYPHS
for cid in {33..255}
do
	C=$(printf "\x$(printf '%x' $cid)")
	CUT="-top 6 -left 8"
	STR="$STR$C";
	WIDTH=`pbmtext "!$C!" | pnmcut -top 6 -height 16 | pnmcrop -left -right | head -n 2 | tail -n 1 | awk '{print $1}' `
	WIDTH=$(echo $WIDTH - 2 | bc)
	PAD=$(echo 16 - $WIDTH | bc)
	pbmtext "!$C!" | pnmcut -top 6 -height 16 | pnmcrop -left -right | pnmcut -left 1 -width $WIDTH | pnmpad -right $PAD -white >../temp/temp.pbm
	convert ../temp/temp.pbm ../temp/temp.xbm
	echo "{ $(grep 0x ../temp/temp.xbm | sed 's/};//') },"  >> $GLYPHS
	#echo "{ $(grep 0x ../temp/temp.xbm) },"  >> $GLYPHS
done
echo "};" >> $GLYPHS

rm -f ../temp/temp.pbm ../temp/temp.xbm

