#include "FlowForce.h"


FlowForce::FlowForce(float v) {
	// Parameters
	float chanLength = 1000.0f;
	float chanWidth = 100.0f;
	float buffArea = 300.0f * 100.0f;
	chanVel0 = v;

	chanHalfLen = 0.5f *chanLength;
	chanHalfWidth = 0.5f * chanWidth;

	// Compute the buffer velocity to have equal flow rates.
	buffVel = 4.0f / 3.0f * chanVel0 * chanWidth * chanWidth / buffArea;
}

Vector3 FlowForce::force(Vector3 r, float diffusion) const {
	if (fabs(r.x) < chanHalfLen) {		// A poiseille flow
		if (fabs(r.y) > chanHalfWidth)
			return Vector3(0.0f);
		float ratio = r.y/chanHalfWidth;
		float vx = chanVel0*(1.0f - ratio*ratio);
		return Vector3(vx / diffusion, 0.0f, 0.0f);
	}
	return Vector3(buffVel/diffusion, 0.0f, 0.0f);
}

