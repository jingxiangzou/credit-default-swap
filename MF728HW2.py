# jxzou@bu.edu
# Jingxiang ZOU
# U93588755
# MSMFT class of 2024
# Yield curve construction
# HW 2
import math
import matplotlib.pyplot as plt
import numpy as np


def fixed_leg(swap_rate, N, fwrs):

    fixed_value = 0
    dntr = 1
    # the start out denominator for t = 0
    # we assume the risk-free rate is continuously compounded
    for n in range(N):
        dntr = dntr * math.exp(-0.5 * fwrs[n])
        fixed_value += 0.5 * swap_rate * dntr

    return fixed_value


def floating_leg(fwrs, N):

    floating_value = 0
    dntr = 1

    for n in range(N):
        dntr = dntr * math.exp(-0.5 * fwrs[n])
        floating_value += 0.5 * fwrs[n] * dntr

    return floating_value


def bisection(swap_rate, N, xsup0, xinf0, epsilon, xlist, plist, xnum):

    d = 1
    x_sup = xsup0
    x_inf = xinf0
    x = 0.1
    wl = []
    for i in range(len(xlist)):
        wl += [xlist[i]] * plist[i]

    while abs(d) > epsilon:

        sr = wl + [x] * xnum
        d = fixed_leg(swap_rate, N, sr) - floating_leg(sr, N)

        if abs(d) < epsilon:
            return x

        if d < 0:
            x_sup = x
            x = x_inf / 2 + x / 2
        else:
            x_inf = x
            x = x_sup / 2 + x / 2

    return 0


if __name__ == "__main__":

    # the numbers that are reused
    xinf0 = 0.0
    xsup0 = 0.1
    epsilon = 0.0000000000001
    pmt_times = [2, 4, 6, 8, 10, 14, 20, 30]
    swap_rates = [2.8438 * 0.01, 3.060 * 0.01, 3.126 * 0.01, 3.144 * 0.01,
                  3.150 * 0.01, 3.169 * 0.01, 3.210 * 0.01, 3.237 * 0.01]

    # question a & b & c
    xy1 = bisection(swap_rates[0], pmt_times[0], xsup0, xinf0, epsilon, [], [], 2)
    xy2 = bisection(swap_rates[1], pmt_times[1], xsup0, xinf0, epsilon, [xy1], [2], 2)
    xy3 = bisection(swap_rates[2], pmt_times[2], xsup0, xinf0, epsilon, [xy1, xy2], [2] * 2, 2)
    xy4 = bisection(swap_rates[3], pmt_times[3], xsup0, xinf0, epsilon, [xy1, xy2, xy3], [2] * 3, 2)
    xy5 = bisection(swap_rates[4], pmt_times[4], xsup0, xinf0, epsilon, [xy1, xy2, xy3, xy4], [2] * 4, 2)
    xy7 = bisection(swap_rates[5], pmt_times[5], xsup0, xinf0, epsilon, [xy1, xy2, xy3, xy4, xy5], [2] * 5, 4)
    xy10 = bisection(swap_rates[6], pmt_times[6], xsup0, xinf0, epsilon, [xy1, xy2, xy3, xy4, xy5, xy7],
                     [2] * 5 + [4], 6)
    xy30 = bisection(swap_rates[7], pmt_times[7], xsup0, xinf0, epsilon, [xy1, xy2, xy3, xy4, xy5, xy7, xy10],
                     [2] * 5 + [4] + [6], 40)

    piecewise1 = [xy1, xy2, xy3, xy4, xy5, xy7, xy10, xy30]
    largef = [(math.exp(x) - 1) for x in piecewise1]

    print("intantaneous forward rates = ", piecewise1)
    print("large forward rates = ", largef)

    # Plot a straight diagonal line with ticked style path
    fig, ax = plt.subplots(figsize=(6, 5))
    plt.title("Swap Rates VS Forward Rates")
    plt.xlabel("tenors")
    plt.ylabel("rates")
    x = ['1y', '2y', '3y', '4y', '5y', '7y', '10y', '30y']
    ax.plot(x, swap_rates, label="swap rates", color="red")
    y = piecewise1
    ax.plot(x, y, label="instantaneous forward rates", color="blue")
    ax.legend()
    plt.show()

    # we proceed to question d
    fwr15y = [xy1] * 2 + [xy2] * 2 + [xy3] * 2 + [xy4] * 2 + [xy5] * 2 + [xy7] * 4 + [xy10] * 6 + [xy30] * 10
    fleg15y = floating_leg(fwr15y, 30)
    swap_rate_15y = fleg15y / fixed_leg(1, 30, fwr15y)
    print("\n\nanswer to question d:")
    print("the 15y swap rate = ", swap_rate_15y)

    # we proceed to question e

    pointnum = np.arange(0.5, 30.5, 0.5)
    pointindex = np.arange(1, 61, 1)
    piecenum = [2, 2, 2, 2, 2, 4, 6, 40]
    fwrcurve = [xy1] * 2 + [xy2] * 2 + [xy3] * 2 + [xy4] * 2 + [xy5] * 2 + [xy7] * 4 + [xy10] * 6 + [xy30] * 40
    zerorates = [0.5 * sum(fwrcurve[:n]) / pointnum[n - 1] for n in pointindex]
    discountf = [math.exp(-1 * zerorates[t - 1] * pointnum[t - 1]) for t in pointindex]
    print("\n\nanswer to question e:\nthe zero rates are:", zerorates)
    print("discount factors:", discountf)
    x1 = np.arange(0.5, 30.5, 0.5)
    plt.plot(x1, zerorates, label="zero rates", color="green")
    x2 = [1, 2, 3, 4, 5, 7, 10, 30]
    plt.plot(x2, swap_rates, label="swap rates", color="red")

    plt.title("Zero rates V.S. Swap Rates")
    plt.xlabel("tenors")
    plt.ylabel("rates")
    plt.legend()
    plt.show()

    # we proceed to question (f)
    # shift all forward rates up 100 bps and recalculate the breakeven swap for each benchmark point
    fwrs = [xy1] * 2 + [xy2] * 2 + [xy3] * 2 + [xy4] * 2 + [xy5] * 2 + [xy7] * 4 + [xy10] * 6 + [xy30] * 40
    tenorsum = [2, 4, 6, 8, 10, 14, 20, 60]
    new_fwrs = [x + 0.01 for x in fwrs]
    # print("new fwrs, ", new_fwrs)
    new_swap_rates = [floating_leg(new_fwrs[:x], x) / fixed_leg(1, x, new_fwrs[:x])
                         for x in tenorsum]
    print("\nanswer to question f:")
    print("the new swap rates = ", new_swap_rates)

    # we plot the new swap rates against the old swap rates
    plt.plot([1, 2, 3, 4, 5, 7, 10, 30], swap_rates, label="original rates", color="red")
    plt.plot([1, 2, 3, 4, 5, 7, 10, 30], new_swap_rates, label="new swap rates", color="purple")

    plt.title("Changes In Swap Rates")
    plt.xlabel("tenors")
    plt.ylabel("rates")
    plt.legend()
    plt.show()

    # g h i j
    bear_swap_rates = [0.028438, 0.03060, 0.03126, 0.03194, 0.03250, 0.03319, 0.03460, 0.03737]
    print(bear_swap_rates)

    brxy1 = bisection(bear_swap_rates[0], pmt_times[0], xsup0, xinf0, epsilon, [], [], 2)
    brxy2 = bisection(bear_swap_rates[1], pmt_times[1], xsup0, xinf0, epsilon, [brxy1], [2], 2)
    brxy3 = bisection(bear_swap_rates[2], pmt_times[2], xsup0, xinf0, epsilon, [brxy1, brxy2], [2] * 2, 2)
    brxy4 = bisection(bear_swap_rates[3], pmt_times[3], xsup0, xinf0, epsilon, [brxy1, brxy2, brxy3], [2] * 3, 2)
    brxy5 = bisection(bear_swap_rates[4], pmt_times[4], xsup0, xinf0, epsilon, [brxy1, brxy2, brxy3, brxy4], [2] * 4, 2)
    brxy7 = bisection(bear_swap_rates[5], pmt_times[5], xsup0, xinf0, epsilon, [brxy1, brxy2, brxy3, brxy4, brxy5], [2] * 5, 4)
    brxy10 = bisection(bear_swap_rates[6], pmt_times[6], xsup0, xinf0, epsilon, [brxy1, brxy2, brxy3, brxy4, brxy5, brxy7],
                     [2] * 5 + [4], 6)
    brxy30 = bisection(bear_swap_rates[7], pmt_times[7], xsup0, xinf0, epsilon, [brxy1, brxy2, brxy3, brxy4, brxy5, brxy7, brxy10],
                     [2] * 5 + [4] + [6], 40)

    piecewise2 = [brxy1, brxy2, brxy3, brxy4, brxy5, brxy7, brxy10, brxy30]
    brlargef = [(math.exp(x) - 1) for x in piecewise2]

    plt.plot(['1y', '2y', '3y', '4y', '5y', '7y', '10y', '30y'], piecewise1, label="original forward rates", color="red")
    plt.plot(['1y', '2y', '3y', '4y', '5y', '7y', '10y', '30y'], piecewise2, label="bear steepened forward rates", color="purple")

    plt.title("Changes In Swap Rates")
    plt.xlabel("tenors")
    plt.ylabel("rates")
    plt.legend()
    plt.show()

    # bull steepener
    bu_swap_rates = [0.023438, 0.02810, 0.02976, 0.03044, 0.03100, 0.03169, 0.03210, 0.03237]
    print(bu_swap_rates)

    buxy1 = bisection(bu_swap_rates[0], pmt_times[0], xsup0, xinf0, epsilon, [], [], 2)
    buxy2 = bisection(bu_swap_rates[1], pmt_times[1], xsup0, xinf0, epsilon, [buxy1], [2], 2)
    buxy3 = bisection(bu_swap_rates[2], pmt_times[2], xsup0, xinf0, epsilon, [buxy1, buxy2], [2] * 2, 2)
    buxy4 = bisection(bu_swap_rates[3], pmt_times[3], xsup0, xinf0, epsilon, [buxy1, buxy2, buxy3], [2] * 3, 2)
    buxy5 = bisection(bu_swap_rates[4], pmt_times[4], xsup0, xinf0, epsilon, [buxy1, buxy2, buxy3, buxy4], [2] * 4, 2)
    buxy7 = bisection(bu_swap_rates[5], pmt_times[5], xsup0, xinf0, epsilon,
                      [buxy1, buxy2, buxy3, buxy4, buxy5], [2] * 5, 4)
    buxy10 = bisection(bu_swap_rates[6], pmt_times[6], xsup0, xinf0, epsilon,
                       [buxy1, buxy2, buxy3, buxy4, buxy5, buxy7], [2] * 5 + [4], 6)
    buxy30 = bisection(bu_swap_rates[7], pmt_times[7], xsup0, xinf0, epsilon, [buxy1, buxy2,
                       buxy3, buxy4, buxy5, buxy7, buxy10], [2] * 5 + [4] + [6], 40)

    piecewise3 = [buxy1, buxy2, buxy3, buxy4, buxy5, buxy7, buxy10, buxy30]
    bulargef = [(math.exp(x) - 1) for x in piecewise3]

    plt.plot(['1y', '2y', '3y', '4y', '5y', '7y', '10y', '30y'], piecewise1, label="original forward rates", color="red")
    plt.plot(['1y', '2y', '3y', '4y', '5y', '7y', '10y', '30y'], piecewise3, label="bull steepened forward rates", color="purple")

    plt.title("Changes In Swap Rates")
    plt.xlabel("tenors")
    plt.ylabel("rates")
    plt.legend()
    plt.show()

































