// Exclude.h
// Copyright Justin Dufresne and Terrance Howard, 2013

#pragma once
#include "useful.h"

struct Restraint {
public:
    Restraint() : id(-1) {}
    Restraint(int id, Vector3 r0, float k) : id(id), r0(r0), k(k) {}
    int id;
    Vector3 r0;
    float k;
};
