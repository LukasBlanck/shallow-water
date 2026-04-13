#pragma once

#include "configs/config.hpp"
#include <filesystem>

SimulationConfig load_config(const std::filesystem::path &path);