/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006 Warren Chou
 Copyright (C) 2007, 2008 StatPro Italia srl

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

#include <ql/instruments/varianceswap.hpp>
#include <ql/event.hpp>

namespace QuantLib {

    VarianceSwap::VarianceSwap(
                          Position::Type position,
                          Real strike,
                          Real notional,
                          const Date& startDate,
                          const Date& maturityDate,
						  const std::vector<Date>& fixingDates,
                          Size daysPerYear,
						  const VarSwapType::Type vswType,
                          const VarSwapBarrierConvention::Type vswConvention,
		                  Real lo_barrier, Real hi_barrier)
    : position_(position), strike_(strike), notional_(notional),
      startDate_(startDate), maturityDate_(maturityDate), 
	  fixingDates_(fixingDates), daysPerYear_(daysPerYear), 
	  vswType_(vswType), vswConvention_(vswConvention), 
	  lo_barrier_(lo_barrier), hi_barrier_(hi_barrier) {}

    Real VarianceSwap::variance() const {
        calculate();
        QL_REQUIRE(variance_ != Null<Real>(), "result not available");
        return variance_;
    }

    void VarianceSwap::setupExpired() const {
        Instrument::setupExpired();
        variance_ = Null<Real>();
    }

    void VarianceSwap::setupArguments(PricingEngine::arguments* args) const {
        VarianceSwap::arguments* arguments =
            dynamic_cast<VarianceSwap::arguments*>(args);
        QL_REQUIRE(arguments != 0, "wrong argument type");

        arguments->position = position_;
        arguments->strike = strike_;
        arguments->notional = notional_;
        arguments->startDate = startDate_;
        arguments->maturityDate = maturityDate_;
        arguments->fixingDates = fixingDates_;
        arguments->daysPerYear = daysPerYear_;
        arguments->vswType = vswType_;
        arguments->vswConvention = vswConvention_;
        arguments->lo_barrier = lo_barrier_;
        arguments->hi_barrier = hi_barrier_;
    }

    void VarianceSwap::fetchResults(const PricingEngine::results* r) const {
        Instrument::fetchResults(r);
        const VarianceSwap::results* results =
            dynamic_cast<const VarianceSwap::results*>(r);
        variance_ = results->variance;
    }

    void VarianceSwap::arguments::validate() const {
        QL_REQUIRE(strike != Null<Real>(), "no strike given");
        QL_REQUIRE(strike > 0.0, "negative or null strike given");
        QL_REQUIRE(notional != Null<Real>(), "no notional given");
        QL_REQUIRE(notional > 0.0, "negative or null notional given");
        QL_REQUIRE(startDate != Date(), "null start date given");
        QL_REQUIRE(maturityDate != Date(), "null maturity date given");
        if (vswType != VarSwapType::Vanilla) {
            QL_REQUIRE(fixingDates.size() > 0,
                       "fixing dates required for non-vanilla variance swap");
            QL_REQUIRE((lo_barrier >=0) 
				|| (hi_barrier>=0), "both lower and upper barriers are invalid");
		}
    }

    bool VarianceSwap::isExpired() const {
        return detail::simple_event(maturityDate_).hasOccurred();
    }

}
