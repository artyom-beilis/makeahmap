#pragma once

class dimage_base {
public:
	static std::string output_directory;
};

class dimage_gray : public dimage_base {
public:
	dimage_gray(std::string const &name,int width,int heigh) :
		file_((output_directory + "/" + name + ".pgm").c_str(),std::fstream::binary);
	{
		file_ << "P5 " << width << " " <<height << " 255\n";
	}
	void write(unsigned char color)
	{
		file_ << color;
	}
	void write(unsigned char r,unsigned char g, unsigned char b)
	{
		file_ << r << g << b;
	}
private:
	std::ofstream file_;
};

class dimage_color : public dimage_base {
public:
	dimage_color(std::string const &name,int width,int height):
		file_((output_directory + "/" + name + ".ppm").c_str(),std::fstream::binary);
	{
		file_ << "P6 " << width << " " <<height << " 255\n";
	}
	void write(unsigned char r,unsigned char g, unsigned char b)
	{
		file_ << r << g << b;
	}
private:
	std::ofstream file_;
};
