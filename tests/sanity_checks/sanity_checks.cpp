#include "tests/sanity_checks/sanity_checks.hpp"

#include "tests/sanity_checks/debug.hpp"
#include "tests/sanity_checks/mass_conservation.hpp"
#include "tests/sanity_checks/positivity.hpp"

#include <memory>
#include <vector>

std::vector<std::unique_ptr<SanityCheck>> make_sanity_checks(const SimulationConfig &cfg) {
    std::vector<std::unique_ptr<SanityCheck>> checks;

    if (cfg.sanity_checks.mass_conservation) {
        checks.push_back(std::make_unique<MassConservation>(cfg.sanity_checks.mass_threshold));
    }
    if (cfg.sanity_checks.positivity) {
        checks.push_back(std::make_unique<Positivity>());
    }
    if (cfg.sanity_checks.debug) {
        checks.push_back(std::make_unique<Debug>());
    }

    return checks;
}