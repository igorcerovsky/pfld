#ifndef POINT_H_
#define POINT_H_


#include <cmath>
#include <limits>

namespace pfld {


template <typename T>
class Point3D
{
public:
	virtual ~Point3D() {};
	Point3D() : x(0), y(0), z(0) { };
	Point3D(T x_, T y_, T z_) : x(x_), y(y_), z(z_) { };
	Point3D(const Point3D& pt) : x(pt.x), y(pt.y), z(pt.z) {}
	const Point3D& operator = (const Point3D& pt)  { x = pt.x; y = pt.y; z = pt.z; return *this; }
	Point3D(const Point3D&& pt)
	{
		x = std::move(pt.x); y = std::move(pt.y); z = std::move(pt.z);
	}
	const Point3D& operator = (const Point3D&& pt)
	{
		x = std::move(pt.x); y = std::move(pt.y); z = std::move(pt.z); return *this;
	}

	bool operator==(const Point3D& pt) const { return (x == pt.x && y == pt.y && z == pt.z); }
	bool operator!=(const Point3D& pt) const { return!(*this == pt); }

	Point3D operator-(const Point3D& pt) const { return Point3D(x - pt.x, y - pt.y, z - pt.z); }
	static void Sub(const Point3D& pt1, const Point3D& pt2, Point3D& res)
	{
		res.x = pt1.x - pt2.x; res.y = pt1.y - pt2.y; res.z = pt1.z - pt2.z;
	}
	Point3D operator+(const Point3D& pt) const { return Point3D(x + pt.x, y + pt.y, z + pt.z); }
	static void Add(const Point3D& pt1, const Point3D& pt2, Point3D& res)
	{
		res.x = pt1.x + pt2.x; res.y = pt1.y + pt2.y; res.z = pt1.z + pt2.z;
	}

	void operator-=(const Point3D& pt) { x -= pt.x; y -= pt.y; z -= pt.z; }
	void operator+=(const Point3D& pt) { x += pt.x; y += pt.y; z += pt.z; }

	// vector arithmetics
	// scalar multiplication of two vectors in 3D
	T operator *(const Point3D& pt) const { return (x*pt.x + y*pt.y + z*pt.z); }
	// vector * number
	Point3D operator *(T a) const { return Point3D(a*x, a*y, a*z); };

	// vector cross product
	Point3D operator /(const Point3D& pt) const
	{
		return Point3D(y*pt.z - z*pt.y,
			z*pt.x - x*pt.z,
			x*pt.y - y*pt.x);
	}
	static void Cross(const Point3D& pt1, const Point3D& pt2, Point3D& res)
	{
		res.x = pt1.y*pt2.z - pt1.z*pt2.y;
		res.y = pt1.z*pt2.x - pt1.x*pt2.z;
		res.z = pt1.x*pt2.y - pt1.y*pt2.x;
	}
	// vector / number
	Point3D operator /(const T a) { return CPoint3D(x / a, y / a, z / a); }

	// vector length
	T Abs() const { return sqrt(x*x + y*y + z*z); };
	// unit vector
	void Unit() {
		T l = Abs();
		if (l > std::numeric_limits<T>().epsilon()) {
			x /= l; y /= l; z /= l;
		}
		else { x = 0; y = 0; z = 0; }
	}


public:
	T x;
	T y;
	T z;
};

} // namespace pfld
#endif