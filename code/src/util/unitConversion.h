#ifndef SBE_WH_UNITCONVERSION_H
#define SBE_WH_UNITCONVERSION_H

namespace atomicUnits{
    // time
    inline double from_zs(const double v) { return v / 2.4188843266049513e4; }
    inline double to_zs(const double v) { return v * 2.4188843266049513e4; }

    inline double from_as(const double v) { return v / 2.4188843266049513e1; }
    inline double to_as(const double v) { return v * 2.4188843266049513e1; }

    inline double from_fs(const double v) { return v / 2.4188843266049513e-2; }
    inline double to_fs(const double v) { return v * 2.4188843266049513e-2; }

    inline double from_ps(const double v) { return v / 2.4188843266049513e-5; }
    inline double to_ps(const double v) { return v * 2.4188843266049513e-5; }

    inline double from_ns(const double v) { return v / 2.4188843266049513e-8; }
    inline double to_ns(const double v) { return v * 2.4188843266049513e-8; }

    inline double from_us(const double v) { return v / 2.4188843266049513e-11; }
    inline double to_us(const double v) { return v * 2.4188843266049513e-11; }

    inline double from_ms(const double v) { return v / 2.4188843266049513e-14; }
    inline double to_ms(const double v) { return v * 2.4188843266049513e-14; }

    inline double from_s(const double v) { return v / 2.4188843266049513e-17; }
    inline double to_s(const double v) { return v * 2.4188843266049513e-17; }

    // length
    inline double from_A(const double v) { return v / 5.291772109060855e-1; }
    inline double to_A(const double v) { return v * 5.291772109060855e-1; }

    inline double from_nm(const double v) { return v / 5.291772109060855e-2; }
    inline double to_nm(const double v) { return v * 5.291772109060855e-2; }

    inline double from_um(const double v) { return v / 5.291772109060855e-5; }
    inline double to_um(const double v) { return v * 5.291772109060855e-5; }

    inline double from_mm(const double v) { return v / 5.291772109060855e-8; }
    inline double to_mm(const double v) { return v * 5.291772109060855e-8; }

    inline double from_m(const double v) { return v / 5.291772109060855e-11; }
    inline double to_m(const double v) { return v * 5.291772109060855e-11; }

    // intensity
    inline double from_W_cm2(const double v) { return v / 3.50944758e16; }
    inline double to_W_cm2(const double v) { return v * 3.50944758e16; }

    inline double from_W_m2(const double v) { return v / 3.50944758e20; }
    inline double to_W_m2(const double v) { return v * 3.50944758e20; }

    // field strength
    inline double from_V_nm(const double v) { return v / 514.22067475617267; }
    inline double to_V_nm(const double v) { return v * 514.22067475617267; }

    inline double from_V_um(const double v) { return v / 514220.67475617267; }
    inline double to_V_um(const double v) { return v * 514220.67475617267; }

    inline double from_V_mm(const double v) { return v / 514220674.75617267; }
    inline double to_V_mm(const double v) { return v * 514220674.75617267; }

    inline double from_V_m(const double v) { return v / 514220674756.17267; }
    inline double to_V_m(const double v) { return v * 514220674756.17267; }

    // energy
    inline double from_eV(const double v) { return v / 27.211386245771678; }
    inline double to_eV(const double v) { return v * 27.211386245771678; }

    inline double from_meV(const double v) { return v / 27211.386245771678; }
    inline double to_meV(const double v) { return v * 27211.386245771678; }

    inline double from_K(const double v) { return v / 315775.0248015561; }
    inline double to_K(const double v) { return v * 315775.0248015561; }

    const double speedOfLight = 137.035999177;
};


#endif
