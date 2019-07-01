/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Warren Chou

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file mcvolatilityhestonengine.hpp
    \brief Monte Carlo variance-swap engine
*/

#ifndef quantlib_mc_volatilityswap_heston_engine_hpp
#define quantlib_mc_volatilityswap_heston_engine_hpp

//#include <ql/pricingengines/vanilla/mcvarianceswapengine.hpp>
#include <ql/instruments/volatilityswap.hpp>
#include <ql/math/integrals/segmentintegral.hpp>
#include <ql/pricingengines/mcsimulation.hpp>
#include <ql/processes/hestonprocess.hpp>
#include <algorithm>
#include <cmath>
//#include <ql/shared_ptr.hpp>

namespace QuantLib {

    //! Variance-swap pricing engine using Monte Carlo simulation,
    /*! as described in Demeterfi, Derman, Kamal & Zou,
        "A Guide to Volatility and Variance Swaps", 1999

        \ingroup forwardengines

        \todo define tolerance of numerical integral and incorporate it
              in errorEstimate

        \test returned fair variances checked for consistency with
              implied volatility curve.
    */
    template <class RNG = PseudoRandom, class S = Statistics>
    class MCVolatilitySwapHestonEngine : public VolatilitySwap::engine,
                                         public McSimulation<MultiVariate, RNG, S> {
      public:
        typedef
            typename McSimulation<MultiVariate, RNG, S>::path_generator_type path_generator_type;
        typedef typename McSimulation<MultiVariate, RNG, S>::path_pricer_type path_pricer_type;
        typedef typename McSimulation<MultiVariate, RNG, S>::stats_type stats_type;
        // constructor
        MCVolatilitySwapHestonEngine(const ext::shared_ptr<HestonProcess>& process,
                                     Size timeSteps,
                                     Size timeStepsPerYear,
                                     bool brownianBridge,
                                     bool antitheticVariate,
                                     Size requiredSamples,
                                     Real requiredTolerance,
                                     Size maxSamples = 999999,
                                     BigNatural seed = 97);
        // calculate variance via Monte Carlo
        void calculate() const {
            McSimulation<MultiVariate, RNG, S>::calculate(requiredTolerance_, requiredSamples_,
                                                          maxSamples_);
            results_.variance = this->mcModel_->sampleAccumulator().mean();

            DiscountFactor riskFreeDiscount =
                process_->riskFreeRate()->discount(arguments_.maturityDate);
            Real multiplier;
            switch (arguments_.position) {
                case Position::Long:
                    multiplier = 1.0;
                    break;
                case Position::Short:
                    multiplier = -1.0;
                    break;
                default:
                    QL_FAIL("Unknown position");
            }
            multiplier *= riskFreeDiscount * arguments_.notional;

            results_.value =
                multiplier * (std::sqrt(std::max(results_.variance, 0.0)) - arguments_.strike);

            if (RNG::allowsErrorEstimate) {
                Real varianceError = this->mcModel_->sampleAccumulator().errorEstimate();
                results_.errorEstimate = multiplier * varianceError;
            }
        }

      protected:
        // McSimulation implementation
        ext::shared_ptr<path_pricer_type> pathPricer() const;
        TimeGrid timeGrid() const;

        ext::shared_ptr<path_generator_type> pathGenerator() const {

            Size dimensions = process_->factors();

            TimeGrid grid = timeGrid();
            typename RNG::rsg_type gen =
                RNG::make_sequence_generator(dimensions * (grid.size() - 1), seed_);

            return ext::shared_ptr<path_generator_type>(
                new path_generator_type(process_, grid, gen, brownianBridge_));
        }
        // data members
        ext::shared_ptr<HestonProcess> process_;
        Size timeSteps_, timeStepsPerYear_;
        Size requiredSamples_, maxSamples_;
        Real requiredTolerance_;
        bool brownianBridge_;
        BigNatural seed_;
    };


    //! Monte Carlo variance-swap engine factory
    template <class RNG = PseudoRandom, class S = Statistics>
    class MakeMCVolatilitySwapHestonEngine {
      public:
        MakeMCVolatilitySwapHestonEngine(const ext::shared_ptr<HestonProcess>& process);
        // named parameters
        MakeMCVolatilitySwapHestonEngine& withSteps(Size steps);
        MakeMCVolatilitySwapHestonEngine& withStepsPerYear(Size steps);
        MakeMCVolatilitySwapHestonEngine& withBrownianBridge(bool b = true);
        MakeMCVolatilitySwapHestonEngine& withSamples(Size samples);
        MakeMCVolatilitySwapHestonEngine& withAbsoluteTolerance(Real tolerance);
        MakeMCVolatilitySwapHestonEngine& withMaxSamples(Size samples);
        MakeMCVolatilitySwapHestonEngine& withSeed(BigNatural seed);
        MakeMCVolatilitySwapHestonEngine& withAntitheticVariate(bool b = true);
        // conversion to pricing engine
        operator ext::shared_ptr<PricingEngine>() const;

      private:
        ext::shared_ptr<HestonProcess> process_;
        bool antithetic_;
        Size steps_, stepsPerYear_, samples_, maxSamples_;
        Real tolerance_;
        bool brownianBridge_;
        BigNatural seed_;
    };

    class VolatilityHestonPathPricer : public PathPricer<MultiPath> {
      public:
        VolatilityHestonPathPricer(const ext::shared_ptr<HestonProcess>& process)
        : process_(process) {}
        Real operator()(const MultiPath& multipath) const;

      private:
        ext::shared_ptr<HestonProcess> process_;
    };

    // inline definitions

    template <class RNG, class S>
    inline MCVolatilitySwapHestonEngine<RNG, S>::MCVolatilitySwapHestonEngine(
        const ext::shared_ptr<HestonProcess>& process,
        Size timeSteps,
        Size timeStepsPerYear,
        bool brownianBridge,
        bool antitheticVariate,
        Size requiredSamples,
        Real requiredTolerance,
        Size maxSamples,
        BigNatural seed)
    : McSimulation<MultiVariate, RNG, S>(antitheticVariate, false), process_(process),
      timeSteps_(timeSteps), timeStepsPerYear_(timeStepsPerYear), requiredSamples_(requiredSamples),
      maxSamples_(maxSamples), requiredTolerance_(requiredTolerance),
      brownianBridge_(brownianBridge), seed_(seed) {
        QL_REQUIRE(timeSteps != Null<Size>() || timeStepsPerYear != Null<Size>(),
                   "no time steps provided");
        QL_REQUIRE(timeSteps == Null<Size>() || timeStepsPerYear == Null<Size>(),
                   "both time steps and time steps per year were provided");
        QL_REQUIRE(timeSteps != 0, "timeSteps must be positive, " << timeSteps << " not allowed");
        QL_REQUIRE(timeStepsPerYear != 0,
                   "timeStepsPerYear must be positive, " << timeStepsPerYear << " not allowed");
    }


    template <class RNG, class S>
    inline TimeGrid MCVolatilitySwapHestonEngine<RNG, S>::timeGrid() const {

        Time t = this->process_->time(this->arguments_.maturityDate);

        if (timeSteps_ != Null<Size>()) {
            return TimeGrid(t, this->timeSteps_);
        } else if (timeStepsPerYear_ != Null<Size>()) {
            Size steps = static_cast<Size>(timeStepsPerYear_ * t);
            return TimeGrid(t, std::max<Size>(steps, 1));
        } else {
            QL_FAIL("time steps not specified");
        }
    }


    template <class RNG, class S>
    inline ext::shared_ptr<typename MCVolatilitySwapHestonEngine<RNG, S>::path_pricer_type>
    MCVolatilitySwapHestonEngine<RNG, S>::pathPricer() const {

        return ext::shared_ptr<typename MCVolatilitySwapHestonEngine<RNG, S>::path_pricer_type>(
            new VolatilityHestonPathPricer(process_));
    }


    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>::MakeMCVolatilitySwapHestonEngine(
        const ext::shared_ptr<HestonProcess>& process)
    : process_(process), antithetic_(false), steps_(Null<Size>()), stepsPerYear_(Null<Size>()),
      samples_(Null<Size>()), maxSamples_(Null<Size>()), tolerance_(Null<Real>()),
      brownianBridge_(false), seed_(0) {}

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withSteps(Size steps) {
        steps_ = steps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withStepsPerYear(Size steps) {
        stepsPerYear_ = steps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withSamples(Size samples) {
        QL_REQUIRE(tolerance_ == Null<Real>(), "tolerance already set");
        samples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withAbsoluteTolerance(Real tolerance) {
        QL_REQUIRE(samples_ == Null<Size>(), "number of samples already set");
        QL_REQUIRE(RNG::allowsErrorEstimate, "chosen random generator policy "
                                             "does not allow an error estimate");
        tolerance_ = tolerance;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withMaxSamples(Size samples) {
        maxSamples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withSeed(BigNatural seed) {
        seed_ = seed;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withBrownianBridge(bool brownianBridge) {
        brownianBridge_ = brownianBridge;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>&
    MakeMCVolatilitySwapHestonEngine<RNG, S>::withAntitheticVariate(bool b) {
        antithetic_ = b;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVolatilitySwapHestonEngine<RNG, S>::
    operator ext::shared_ptr<PricingEngine>() const {
        QL_REQUIRE(steps_ != Null<Size>() || stepsPerYear_ != Null<Size>(),
                   "number of steps not given");
        QL_REQUIRE(steps_ == Null<Size>() || stepsPerYear_ == Null<Size>(),
                   "number of steps overspecified");
        return ext::shared_ptr<PricingEngine>(new MCVolatilitySwapHestonEngine<RNG, S>(
            process_, steps_, stepsPerYear_, brownianBridge_, antithetic_, samples_, tolerance_,
            maxSamples_, seed_));
    }


    inline Real VolatilityHestonPathPricer::operator()(const MultiPath& multiPath) const {
        // const Path& path = multiPath[0];
        QL_REQUIRE(multiPath.pathSize() > 0, "the path cannot be empty");
        Time t0 = multiPath[0].timeGrid().front();
        Time t = multiPath[0].timeGrid().back();
        Time dt = multiPath[0].timeGrid().dt(0);
        SegmentIntegral integrator(static_cast<Size>(t / dt));
        detail_heston::Integrand f(multiPath, process_);
        return integrator(f, t0, t) / t;
    }
}


#endif
