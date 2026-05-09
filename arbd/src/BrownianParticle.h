// BrownianParticle.h
// Copyright Justin Dufresne and Terrance Howard, 2013

#include <algorithm>

#include <cuda.h>

#include "useful.h"
#include "JamesBond.h"

class BrownianParticle {
public:
	BrownianParticle() :
			id(-1), type(-1),
			has_orientation_(false) { }

	BrownianParticle(int id) :
			id(id), type(0), pos(Vector3(0.0f)),
			has_orientation_(false) { }

	BrownianParticle(int id, const Vector3& pos, int type) :
			id(id), type(type), pos(pos),
			has_orientation_(false) { }

	BrownianParticle(int id, const Vector3& pos, int type,
									 const Vector3& orientation) :
			id(id), type(type), pos(pos),
			orientation(orientation),
			has_orientation_(true) { }

	HOST DEVICE
	inline bool has_orientation() const { return has_orientation_; }

	HOST DEVICE
	inline bool is_dummy() const { return is_dummy_; }
	
	HOST DEVICE
	inline void lose_orientation() { has_orientation_ = false; }

	HOST DEVICE
	inline void add_orientation(const Vector3& o) {
		orientation = o;
		has_orientation_ = true;
	}

	BrownianParticle& operator=(const BrownianParticle& src);

	// Static comparison functions for sorting
	static inline bool compareByIndex(const BrownianParticle& a,
																		const BrownianParticle& b) {
		return a.id < b.id;
	}

	static inline bool compareByType(const BrownianParticle& a,
																	 const BrownianParticle& b) {
		if (a.type == b.type)
			return compareByIndex(a, b);
		return a.type < b.type;
	}

public:
	int id;
	int type; // index into global type array.

	Vector3 pos;

	Vector3 orientation;

private:
	bool is_dummy_;

	bool has_orientation_;
};
