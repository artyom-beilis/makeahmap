#pragma once
#include "glyphs.h"
#include <assert.h>

static const int font_digit_width  = FONT_WIDTH;
static const int font_digit_height = FONT_HEIGHT;

template<typename T>
struct get_size;

template<typename M>
struct get_size<std::vector<std::vector<M> > > {
	static int height(std::vector<std::vector<M> > const &m)
	{
		return m.size();
	}
	static int width(std::vector<std::vector<M> > const &m)
	{
		return m.at(0).size();
	}
};

template<typename T>
class image;

template<typename M>
struct get_size<image<M> > {
	static int height(image<M> const &m)
	{
		return m.height();
	}
	static int width(image<M> const &m)
	{
		return m.width();
	}
};


template<typename ImgType>
int print_char(char ch,int r,int c,int color,ImgType *img)
{
	if(img) {
		if(r < 0)
			return 0;
		if(c < 0)
			return 0;
		if(r + font_digit_height >= get_size<ImgType>::height(*img))
			return 0;
		if(c + font_digit_width  >= get_size<ImgType>::width(*img))
			return 0;
	}
	unsigned char uc=ch;
	int glyphs_total = sizeof(font_digits) / sizeof(font_digits[0]);
	int id = int(uc) - 33;
	if(id < 0 || id >= glyphs_total)
		return 0;
	unsigned char *digit = font_digits[id];
	int cols = 0;
	int x_fact = (font_digit_width + 7)/8;
	for(int dr=0;dr<font_digit_height;dr++) {
		for(int dc=0;dc<font_digit_width;dc++) {
			if(digit[dr * x_fact + dc / 8] & (1<<(dc % 8))) {
				if(img) {
					(*img)[r+dr][c+dc]=color;
				}
				cols = std::max(cols,dc + 1);
			}
		}
	}
	return cols;

}

template<typename ImgType>
int print_str_internal(char const *msg,int r,int c,int color,ImgType *img)
{
	int col_start = c;
	while(*msg!=0) {
		if(*msg == ' ') {
			c+=4;
		}
		else {
			int cols = print_char(*msg,r,c,color,img);
			c+=cols + 1;
		}
		msg++;
	}
	return c-col_start;
}

template<typename ImgType>
void print_str(char const *msg,int r,int c,int color,ImgType &img)
{
	print_str_internal<ImgType>(msg,r,c,color,&img);
}

int get_print_str_len(char const *msg)
{
	return print_str_internal<std::vector<std::vector<int> > >(msg,0,0,1,0);
}

