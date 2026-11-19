import math


class BioCalculations:
    @staticmethod
    def egfr_ckd_epi_2021(creatinine_mg_dl: float, age: int, female: bool) -> float:
        if not math.isfinite(creatinine_mg_dl):
            raise ValueError("creatinine must be finite")
        if creatinine_mg_dl <= 0:
            raise ValueError("creatinine must be positive and non-zero")
        if isinstance(age, bool) or not isinstance(age, int):
            raise ValueError("age must be an integer")
        if age < 18 or age > 120:
            raise ValueError("Adult CKD-EPI equation is validated for ages 18 to 120")
        kappa = 0.7 if female else 0.9
        alpha = -0.241 if female else -0.302
        sex_factor = 1.012 if female else 1.0
        ratio = creatinine_mg_dl / kappa
        return 142 * (min(ratio, 1) ** alpha) * (max(ratio, 1) ** -1.200) * (0.9938 ** age) * sex_factor

    @staticmethod
    def creatinine_cockcroft_gault(creatinine_mg_dl: float, age: int, weight_kg: float, female: bool) -> float:
        if not math.isfinite(creatinine_mg_dl) or not math.isfinite(weight_kg):
            raise ValueError("creatinine and weight must be finite")
        if creatinine_mg_dl <= 0 or weight_kg <= 0:
            raise ValueError("creatinine and weight must be positive")
        if weight_kg > 500:
            raise ValueError("weight exceeds supported adult range")
        if isinstance(age, bool) or not isinstance(age, int):
            raise ValueError("age must be an integer")
        if age < 18 or age > 120:
            raise ValueError("Cockcroft-Gault formula is validated for adult ages 18 to 120")
        result = ((140 - age) * weight_kg) / (72 * creatinine_mg_dl)
        return result * 0.85 if female else result

    @staticmethod
    def anion_gap(sodium: float, chloride: float, bicarbonate: float, albumin: float | None = None) -> float:
        values = (sodium, chloride, bicarbonate)
        if any(not math.isfinite(v) for v in values):
            raise ValueError("electrolyte values must be finite")
        gap = sodium - (chloride + bicarbonate)
        if albumin is not None:
            if not math.isfinite(albumin):
                raise ValueError("albumin must be finite")
            gap += 2.5 * (4.0 - albumin)
        return gap

    @staticmethod
    def hepatocyte_r_ratio(alt: float, alt_uln: float, alp: float, alp_uln: float) -> float:
        if any(not math.isfinite(value) for value in (alt, alt_uln, alp, alp_uln)):
            raise ValueError("R-ratio values must be finite")
        if alt < 0 or alt_uln <= 0 or alp <= 0 or alp_uln <= 0:
            raise ValueError("ALT/ALP must be non-negative, and ALP/ULN values must be strictly positive non-zero")
        return (alt / alt_uln) / (alp / alp_uln)

    @staticmethod
    def ifcc_to_ngsp(ifcc_mmol_mol: float) -> float:
        if not math.isfinite(ifcc_mmol_mol):
            raise ValueError("HbA1c must be finite")
        return 0.09148 * ifcc_mmol_mol + 2.152
