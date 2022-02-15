def get_freq(fname: str, unit: str = "MHz") -> float:
    return float(fname.split("Freq")[-1].split(f"{unit}")[0])
