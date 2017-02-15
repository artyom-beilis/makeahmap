#pragma once
#include <stdlib.h>
template<typename T>
class image {
public:
	image(int w = 0,int h = 0,int wpad = 0,int hpad = 0) :
		width_(w),
		height_(h),
		wpad_(wpad),
		hpad_(hpad),
		data_(nullptr)
	{
		mem_width_ = wpad_*2 +  width_;
		mem_size_  = mem_width_ * (h+hpad*2);
		data_ = reinterpret_cast<T*>(allocate(mem_size_ * sizeof(T)));
	}
	image(image const &)=delete;
	void operator=(image const &) = delete;
	
	image (image &&o) 
	{
		move_and_reset(o);
	}

	

	image &operator=(image &&o) 
	{
		if(this != &o) {
			move_and_reset(o);
		}
		return *this;
	}
	
	int height() const { return height_; }
	int width() const { return width_; }

	T const *operator[](int index) const
	{
		return data_ + (index + hpad_)*mem_width_ + wpad_;
	}
	T *operator[](int index)
	{
		return data_ + (index + hpad_)*mem_width_ + wpad_;
	}

	~image()
	{
		deallocate(data_);
	}
	
private:
	void move_and_reset(image &o)
	{
		width_ = o.width_;
		height_ = o.height_;
		wpad_ = o.wpad_;
		hpad_ = o.hpad_;
		mem_width_ = o.mem_width_;
		mem_size_  = o.mem_size_;
		data_ = o.data_;
		memset(&o,0,sizeof(o));
	}
	static void *allocate(size_t size)
	{
		if(size == 0)
			return nullptr;
		#if defined _WIN32 || defined WIN32
		void *ptr = _aligned_malloc(size,32);
		#else	
		void *ptr = aligned_alloc(32,size);
		#endif
		if(!ptr)
			throw std::bad_alloc();
		memset(ptr,0,size);
		return ptr;
	}
	static void deallocate(void *p)
	{
		#if defined _WIN32 || defined WIN32
		if(p)
			_aligned_free(p);
		#else
		free(p);
		#endif
	}
	int width_;
	int height_;
	int wpad_;
	int hpad_;
	int mem_width_;
	int mem_size_;
	T *data_;
};
