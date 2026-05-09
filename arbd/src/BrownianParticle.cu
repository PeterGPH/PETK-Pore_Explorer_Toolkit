// BrownianParticle.cu
// Copyright Justin Dufresne and Terrance Howard, 2013

#include "BrownianParticle.h"

BrownianParticle& BrownianParticle::operator=(const BrownianParticle& src) {
	id = src.id;
	type = src.type;
	pos = src.pos;
	orientation = src.orientation;
	has_orientation_ = src.has_orientation_;
	is_dummy_ = src.is_dummy_;
	return *this;
}
