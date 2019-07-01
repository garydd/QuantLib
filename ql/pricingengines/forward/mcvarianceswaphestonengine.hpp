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

/*! \file mcvarianceswaphestonengine.hpp
    \brief Monte Carlo variance-swap engine
*/

#ifndef quantlib_mc_varianceswap_heston_engine_hpp
#define quantlib_mc_varianceswap_heston_engine_hpp

//#include <ql/pricingengines/vanilla/mcvarianceswapengine.hpp>
#include <ql/instruments/varianceswap.hpp>
#include <ql/math/integrals/segmentintegral.hpp>
#include <ql/pricingengines/mcsimulation.hpp>
#include <ql/processes/hestonprocess.hpp>
//#include <ql/shared_ptr.hpp>
#include <iostream>

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
    class MCVarianceSwapHestonEngine : public VarianceSwap::engine,
                                       public McSimulation<MultiVariate, RNG, S> {
      public:
        typedef
            typename McSimulation<MultiVariate, RNG, S>::path_generator_type path_generator_type;
        typedef typename McSimulation<MultiVariate, RNG, S>::path_pricer_type path_pricer_type;
        typedef typename McSimulation<MultiVariate, RNG, S>::stats_type stats_type;
        // constructor
        MCVarianceSwapHestonEngine(const ext::shared_ptr<HestonProcess>& process,
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

            if (this->arguments_.fixingDates.size() > 0) {
                results_.variance *= this->arguments_.daysPerYear;
            }

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

            results_.value = multiplier * (results_.variance - arguments_.strike);

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
    class MakeMCVarianceSwapHestonEngine {
      public:
        MakeMCVarianceSwapHestonEngine(const ext::shared_ptr<HestonProcess>& process);
        // named parameters
        MakeMCVarianceSwapHestonEngine& withSteps(Size steps);
        MakeMCVarianceSwapHestonEngine& withStepsPerYear(Size steps);
        MakeMCVarianceSwapHestonEngine& withBrownianBridge(bool b = true);
        MakeMCVarianceSwapHestonEngine& withSamples(Size samples);
        MakeMCVarianceSwapHestonEngine& withAbsoluteTolerance(Real tolerance);
        MakeMCVarianceSwapHestonEngine& withMaxSamples(Size samples);
        MakeMCVarianceSwapHestonEngine& withSeed(BigNatural seed);
        MakeMCVarianceSwapHestonEngine& withAntitheticVariate(bool b = true);
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

    class VarianceHestonPathPricer : public PathPricer<MultiPath> {
      public:
        VarianceHestonPathPricer(const ext::shared_ptr<HestonProcess>& process)
        : process_(process) {}
        Real operator()(const MultiPath& multipath) const;

      private:
        ext::shared_ptr<HestonProcess> process_;
    };

    class VarianceBarrierHestonPathPricer : public PathPricer<MultiPath> {
      public:
        VarianceBarrierHestonPathPricer(const ext::shared_ptr<HestonProcess>& process,
                                        const VarSwapType::Type vswType,
                                        const VarSwapBarrierConvention::Type vswConvention,
                                        const Real lo_barrier,
                                        const Real hi_barrier)
        : process_(process), vswType_(vswType), vswConvention_(vswConvention), 
			lo_barrier_(lo_barrier), hi_barrier_(hi_barrier) {}
        Real operator()(const MultiPath& multipath) const;

      private:
        ext::shared_ptr<HestonProcess> process_;
        const VarSwapType::Type vswType_;
        const VarSwapBarrierConvention::Type vswConvention_;
        const Real lo_barrier_;
        const Real hi_barrier_;
    };
    // inline definitions

    template <class RNG, class S>
    inline MCVarianceSwapHestonEngine<RNG, S>::MCVarianceSwapHestonEngine(
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
    inline TimeGrid MCVarianceSwapHestonEngine<RNG, S>::timeGrid() const {
        Time t = this->process_->time(this->arguments_.maturityDate);
        std::vector<Time> fixing_times;
        // QL_REQUIRE(!(this->arguments_.fixingDates.empty()), "empty fixing dates")
        for (std::vector<Date>::iterator f_date = this->arguments_.fixingDates.begin();
             f_date != this->arguments_.fixingDates.end(); ++f_date) {
            fixing_times.push_back(this->process_->time(*f_date));
            // std::cout << "timegrid" << *f_date << ' ' << this->process_->time(*f_date) << ' \n';
        }
        if (!fixing_times.empty()) {
            QL_REQUIRE(this->arguments_.fixingDates.back() == this->arguments_.maturityDate,
                       "last fixing date" << this->arguments_.fixingDates.back()
                                          << " must be equal to maturity "
                                          << this->arguments_.maturityDate);
        }

        if (timeSteps_ != Null<Size>()) {
            if (!fixing_times.empty()) {
                return TimeGrid(fixing_times.begin(), fixing_times.end(), this->timeSteps_);
            } else {
                return TimeGrid(t, this->timeSteps_);
            }
        } else if (timeStepsPerYear_ != Null<Size>()) {
            Size steps = static_cast<Size>(timeStepsPerYear_ * t);
            if (!fixing_times.empty()) {
                return TimeGrid(fixing_times.begin(), fixing_times.end(), std::max<Size>(steps, 1));
            } else {
                return TimeGrid(t, std::max<Size>(steps, 1));
            }
        } else {
            QL_FAIL("time steps not specified");
        }
    }


    template <class RNG, class S>
    inline ext::shared_ptr<typename MCVarianceSwapHestonEngine<RNG, S>::path_pricer_type>
    MCVarianceSwapHestonEngine<RNG, S>::pathPricer() const {
        if (this->arguments_.vswType == VarSwapType::Vanilla) {
            return ext::shared_ptr<typename MCVarianceSwapHestonEngine<RNG, S>::path_pricer_type>(
                new VarianceHestonPathPricer(process_));
        } else {
            return ext::shared_ptr<typename MCVarianceSwapHestonEngine<RNG, S>::path_pricer_type>(
                new VarianceBarrierHestonPathPricer(process_, 
													this->arguments_.vswType,
													this->arguments_.vswConvention,
													this->arguments_.lo_barrier,
													this->arguments_.hi_barrier));
        }
    }


    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>::MakeMCVarianceSwapHestonEngine(
        const ext::shared_ptr<HestonProcess>& process)
    : process_(process), antithetic_(false), steps_(Null<Size>()), stepsPerYear_(Null<Size>()),
      samples_(Null<Size>()), maxSamples_(Null<Size>()), tolerance_(Null<Real>()),
      brownianBridge_(false), seed_(0) {}

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withSteps(Size steps) {
        steps_ = steps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withStepsPerYear(Size steps) {
        stepsPerYear_ = steps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withSamples(Size samples) {
        QL_REQUIRE(tolerance_ == Null<Real>(), "tolerance already set");
        samples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withAbsoluteTolerance(Real tolerance) {
        QL_REQUIRE(samples_ == Null<Size>(), "number of samples already set");
        QL_REQUIRE(RNG::allowsErrorEstimate, "chosen random generator policy "
                                             "does not allow an error estimate");
        tolerance_ = tolerance;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withMaxSamples(Size samples) {
        maxSamples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withSeed(BigNatural seed) {
        seed_ = seed;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withBrownianBridge(bool brownianBridge) {
        brownianBridge_ = brownianBridge;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>&
    MakeMCVarianceSwapHestonEngine<RNG, S>::withAntitheticVariate(bool b) {
        antithetic_ = b;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCVarianceSwapHestonEngine<RNG, S>::operator ext::shared_ptr<PricingEngine>() const {
        QL_REQUIRE(steps_ != Null<Size>() || stepsPerYear_ != Null<Size>(),
                   "number of steps not given");
        QL_REQUIRE(steps_ == Null<Size>() || stepsPerYear_ == Null<Size>(),
                   "number of steps overspecified");
        return ext::shared_ptr<PricingEngine>(new MCVarianceSwapHestonEngine<RNG, S>(
            process_, steps_, stepsPerYear_, brownianBridge_, antithetic_, samples_, tolerance_,
            maxSamples_, seed_));
    }


    namespace detail_heston {

        class Integrand {
          public:
            Integrand(const MultiPath& multiPath, const ext::shared_ptr<HestonProcess>& process)
            : path_(multiPath), process_(process) {}
            Real operator()(Time t) const {
                Size i = static_cast<Size>(t / path_[0].timeGrid().dt(0));
                /*Array x(path_.pathSize());
                for (int j = 0; j++; j < path_.pathSize()) {
                    x[j] = path_[j][i];
                                }
                Real sigma = process_->diffusion(t, x)[0][0];
                return sigma*sigma;*/
                return path_[1][i]; //*path_[1][i];
            }

          private:
            MultiPath path_;
            ext::shared_ptr<HestonProcess> process_;
        };
    }


    // inline Real VarianceHestonPathPricer::operator()(const MultiPath& multiPath) const {
    //    // const Path& path = multiPath[0];
    //    QL_REQUIRE(multiPath.pathSize() > 0, "the path cannot be empty");
    //    Time t0 = multiPath[0].timeGrid().front();
    //    Time t = multiPath[0].timeGrid().back();
    //    Time dt = multiPath[0].timeGrid().dt(0);
    //    SegmentIntegral integrator(static_cast<Size>(t / dt));
    //    detail_heston::Integrand f(multiPath, process_);
    //    return integrator(f, t0, t) / t;
    //}

    inline Real VarianceHestonPathPricer::operator()(const MultiPath& multiPath) const {
        // const Path& path = multiPath[0];
        QL_REQUIRE(multiPath.pathSize() > 0, "the path cannot be empty");
        Real rv = 0.0;
        TimeGrid tg = multiPath[0].timeGrid();
        std::vector<Time> mtg = tg.mandatoryTimes();
        Size idx, idx_prev;
        Real logReturn;
        int N = 0;
        if (mtg.size() > 1) {
            for (std::vector<Time>::iterator t = mtg.begin(); t != mtg.end(); ++t) {
                idx = tg.closestIndex(*t);
                // std::cout << *t << ' ' << idx << ' \n' ;
                // QL_REQUIRE(multiPath[0].value(idx) == 0.0, "abc" << *(mtg.begin()) << "xyz" <<
                // *(mtg.end()) << "def") <<std::endl;
                if (idx > 0) {
                    logReturn = std::logf(multiPath[0].value(idx) / multiPath[0].value(idx_prev));
                    rv += logReturn * logReturn;
                    N += 1;
                }
                idx_prev = idx;
            }
            // return logReturn;
            // std::cout << "RV: " << rv<<" size" <<mtg.size() <<std::endl;
            return rv / N;
            //(mtg.size() - 1);
        } else {
            // only one possibility here due to timegrid construction above
            Time t0 = multiPath[0].timeGrid().front();
            Time t = multiPath[0].timeGrid().back();
            Time dt = multiPath[0].timeGrid().dt(0);
            SegmentIntegral integrator(static_cast<Size>(t / dt));
            detail_heston::Integrand f(multiPath, process_);
            // std::cout << "t0: " << t0 << " t:" << t << " dt:" << dt
            //<< " size:" << multiPath[0].timeGrid().size() <<std::endl;
            return integrator(f, t0, t) / t;
        }
    }

    inline Real VarianceBarrierHestonPathPricer::operator()(const MultiPath& multiPath) const {
        QL_REQUIRE(multiPath.pathSize() > 0, "the path cannot be empty");
        Real rv = 0.0;
        TimeGrid tg = multiPath[0].timeGrid();
        std::vector<Time> mtg = tg.mandatoryTimes();
        Size idx, idx_prev;
        Real logReturn;
        int N = 0;
        for (std::vector<Time>::iterator t = mtg.begin(); t != mtg.end(); ++t) {
            idx = tg.closestIndex(*t);
            //std::cout << " " << multiPath[0].value(idx) << std::endl;
            if (idx > 0) {
                logReturn = std::logf(multiPath[0].value(idx) / multiPath[0].value(idx_prev));
                bool breach_lo = (this->lo_barrier_ >= 0) ? 
					(multiPath[0].value(idx) < this->lo_barrier_) : false;
                bool breach_hi = (this->hi_barrier_ >= 0) ?
                    (multiPath[0].value(idx) > this->hi_barrier_) : false;
                bool breach_lo_prev = (this->lo_barrier_ >= 0) ?
                    (multiPath[0].value(idx_prev) < this->lo_barrier_) : false;
                bool breach_hi_prev = (this->hi_barrier_ >= 0) ?
                    (multiPath[0].value(idx_prev) > this->hi_barrier_) : false;
                bool breach_overall;
				if (this->vswConvention_ == VarSwapBarrierConvention::N) {
                    breach_overall = breach_lo || breach_hi;
				} else if (this->vswConvention_ == VarSwapBarrierConvention::N_1) {
                    breach_overall = breach_lo_prev || breach_hi_prev;
				} else {
                    breach_overall = breach_lo_prev || breach_hi_prev || breach_lo || breach_hi;
                }
                
				if (!breach_overall) {
                    rv += logReturn * logReturn;
                    N += 1;
                } else {
                    if (this->vswType_ == VarSwapType::Corridor) {
                        N += 1;
					}        
				}
			}
            idx_prev = idx;
        }
        //std::cout << " " << rv <<" "<< N << std::endl;
        return (rv / N) ;
	}
}


#endif
