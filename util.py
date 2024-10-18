def in_bin_freq(bins, vals, probs):
    freqs = [0 for i in range(len(bins))]
    vals.sort()
    for val, prob in zip(vals, probs):
        for i, bin in enumerate(bins):
            if val <= bin:
                freqs[i] += prob
                break

    return freqs