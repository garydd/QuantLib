/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Warren Chou
 Copyright (C) 2008 StatPro Italia srl

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

/*! \file varianceswap.hpp
    \brief Variance swap
*/

#ifndef quantlib_variance_swap_hpp
#define quantlib_variance_swap_hpp

#include <ql/instruments/payoffs.hpp>
#include <ql/option.hpp>
#include <ql/position.hpp>
#include <ql/processes/blackscholesprocess.hpp>

namespace QuantLib {

    //! Variance swap
    /*! \warning This class does not manage seasoned variance swaps.

        \ingroup instruments
    */
    struct VarSwapType {
        enum Type { Vanilla = 0, Corridor, Conditional };
    };

    struct VarSwapBarrierConvention {
        enum Type { N = 0, N_1, N_1_N };
    };

    class VarianceSwap : public Instrument {
      public:
        class arguments;
        class results;
        class engine;

        VarianceSwap(Position::Type position,
                     Real strike,
                     Real notional,
                     const Date& startDate,
                     const Date& maturityDate,
                     const std::vector<Date>& fixingDates,
                     Size daysPerYear,
                     const VarSwapType::Type vswType,
                     const VarSwapBarrierConvention::Type vswConvention,
                     const Real lo_barrier = -1,
					 const Real hi_barrier = -1);
        //! \name Instrument interface
        //@{
        bool isExpired() const;
        //@}
        //! \name Additional interface
        //@{
        // inspectors
        Real strike() const;
        Position::Type position() const;
        Date startDate() const;
        Date maturityDate() const;
        std::vector<Date> fixingDates() const;
        Real notional() const;
        // results
        Real variance() const;

        Size daysPerYear() const;
        VarSwapType::Type vswType() const;
        VarSwapBarrierConvention::Type vswConvention() const;
        Real lo_barrier() const;
        Real hi_barrier() const;
        //@}
        // other
        void setupArguments(PricingEngine::arguments* args) const;
        void fetchResults(const PricingEngine::results*) const;

      protected:
        void setupExpired() const;
        // data members
        Position::Type position_;
        Real strike_;
        Real notional_;
        Date startDate_, maturityDate_;
        std::vector<Date> fixingDates_;
        // results
        mutable Real variance_;
        Size daysPerYear_;
        VarSwapType::Type vswType_;
        VarSwapBarrierConvention::Type vswConvention_;
        Real lo_barrier_;
        Real hi_barrier_;
    };


    //! %Arguments for forward fair-variance calculation
    class VarianceSwap::arguments : public virtual PricingEngine::arguments {
      public:
        arguments() : strike(Null<Real>()), notional(Null<Real>()) {}
        void validate() const;
        Position::Type position;
        Real strike;
        Real notional;
        Date startDate;
        Date maturityDate;
        std::vector<Date> fixingDates;
        Size daysPerYear;
        VarSwapType::Type vswType;
        VarSwapBarrierConvention::Type vswConvention;
        Real lo_barrier;
        Real hi_barrier;
    };


    //! %Results from variance-swap calculation
    class VarianceSwap::results : public Instrument::results {
      public:
        Real variance;
        void reset() {
            Instrument::results::reset();
            variance = Null<Real>();
        }
    };

    //! base class for variance-swap engines
    class VarianceSwap::engine
    : public GenericEngine<VarianceSwap::arguments, VarianceSwap::results> {};


    // inline definitions

    inline Date VarianceSwap::startDate() const { return startDate_; }

    inline Date VarianceSwap::maturityDate() const { return maturityDate_; }

    inline std::vector<Date> VarianceSwap::fixingDates() const { return fixingDates_; }

    inline Real VarianceSwap::strike() const { return strike_; }

    inline Real VarianceSwap::notional() const { return notional_; }

    inline Position::Type VarianceSwap::position() const { return position_; }

    inline Size VarianceSwap::daysPerYear() const { return daysPerYear_; }

    inline VarSwapType::Type VarianceSwap::vswType() const { return vswType_; }

	inline VarSwapBarrierConvention::Type VarianceSwap::vswConvention() const { return vswConvention_; }

	inline Real VarianceSwap::lo_barrier() const {
            return lo_barrier_;
        }

	inline Real VarianceSwap::hi_barrier() const { return hi_barrier_; }
}


#endif
