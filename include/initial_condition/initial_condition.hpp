#pragma once

class Grid; // here only references are needed
class State;

class InitialCondition {
  public:
    virtual ~InitialCondition() = default;
    virtual void apply(const Grid &grid, State &U) const = 0;
};