import math
import matplotlib.pyplot as plt


# we assume that the hazard rate is piece-wise constant between these instances: 0y, 1y, 2y, 3y, 5y


def premium_leg(spread, N, dntrs, intervals, hazard_rates):
    """
    :param spread: corresponding spread, constant for one specific CDS contract.
    :param dntrs: corresponding denominator, list, length = N.
    :param N: the number of premium payments, constant for one specific CDS contract.
    :param intervals: corresponding interval, list, length = N.
    :param hazard_rates : hazard rates, list, length is piece_wise constant, piece-wise constant
    :return: the value of the premium leg
    """
    premiums = 0
    # this in the end returns the value of the premium legs

    prob_begin = 1.0
    # at inception, the survival probability is 1.
    # prb_begin: the probability of not default till the beginning of a specific interval
    # pro_end: the probability of not default till the beginning of a specific interval

    for n in range(1, N + 1, 1):
        prob_tail = prob_begin * math.exp(-1 * hazard_rates[n - 1] * intervals[n - 1])
        premiums += 0.5 * spread * pow(dntrs[n - 1], -1 * n) * intervals[n - 1] * (prob_begin + prob_tail)
        prob_begin = prob_tail

    # the length of hazard rates should be the same as that of the interval

    return premiums


def protection_leg(recovery_rate, N, dntrs, intervals, hazard_rates):
    """
    :param recovery_rate: recovery rate
    :param dntrs: corresponding denominator, list, length = N
    :param N: the number of premium payments, constant for one specific CDS contract.
    :param intervals: corresponding interval, list, length = N
    :param hazard_rates : hazard rates, length = N, piece-wise constant
    :return: the value of the protection leg
    """

    protection = 0
    prob_begin = 1.0

    for n in range(1, N + 1, 1):
        prob_tail = prob_begin * math.exp(-1 * hazard_rates[n - 1] * intervals[n - 1])
        protection += (1 - recovery_rate) * pow(dntrs[n - 1], -1 * n) * (prob_begin - prob_tail)
        prob_begin = prob_tail

    return protection


if __name__ == "__main__":

    # bootstrapping
    print("# bootstrapping step 1", "\n")
    # assume premiums are paid semi-annually
    # for the maturity = 1 contract
    N_1y = 2
    intervals_1y = [0.5, 0.5]
    dntr_1y = [1.03, 1.03]
    d = 1
    r = 0.4  # recovery rate
    x_inf = 0.0  # hazard rate trials starting points
    x_sup = 1.0
    x1 = 0.25
    spread_1y = 0.01

    # bisection method to compute hazard rate of stage 1
    while abs(d) > 0.0000000001:

        hr = [x1, x1]
        # print("x1=", x1)
        d = protection_leg(r, N_1y, dntr_1y, intervals_1y, hr) - \
            premium_leg(spread_1y, N_1y, dntr_1y, intervals_1y, hr)
        # print(d)
        if d > 0:
            x_sup = x1
            x1 = x_inf * 0.5 + x1 * 0.5
            # if the protection leg is over-priced, it means that the hazard rates are over-priced
        else:
            x_inf = x1
            x1 = x_sup * 0.5 + x1 * 0.5

            # if the protection leg is under-priced, it means that the hazard rates are under-priced

    # maturity = 1y CDS:
    print("1y CDS with nominal = $1 : premium leg = ", premium_leg(spread_1y, N_1y, dntr_1y, intervals_1y, [x1, x1]))
    print("1y CDS with nominal = $1 : protection leg = ", protection_leg(r, N_1y, dntr_1y, intervals_1y, [x1, x1]))
    print("hazard rate of 0y ~ 1y = ", x1)
    Prob_0point5y = 1 * math.exp(-1 * x1 * 0.5)
    Prob_1y = 1 * math.exp(-1 * x1 * 1)
    print("prob survival till 0.5y = ", Prob_0point5y)
    print("prob survival till 1y = ", Prob_1y)

    # Bootstrapping step 2
    print("\n\n# bootstrapping step 2", "\n")
    # assume premiums are paid semi-annually
    # for the maturity = 2 contract
    N_2y = 4
    intervals_2y = [0.5, 0.5, 0.5, 0.5]
    dntr_2y = [1.03, 1.03, 1.03, 1.03]
    d = 1
    r = 0.4  # recovery rate
    x_inf = 0.0  # hazard rate trials starting points
    x_sup = 1.0
    x2 = 0.25
    spread_2y = 0.011
    d = 1

    # bisection method to compute hazard rate of stage 2
    while abs(d) > 0.00000001:

        hr = [x1, x1, x2, x2]
        d = protection_leg(r, N_2y, dntr_2y, intervals_2y, hr) - \
            premium_leg(spread_2y, N_2y, dntr_2y, intervals_2y, hr)

        if abs(d) < 0.00000001:
            break

        if d > 0:
            x_sup = x2
            x2 = x_inf * 0.5 + x2 * 0.5
            # if the protection leg is over-priced, it means that the hazard rates are over-priced
        else:
            x_inf = x2
            x2 = x_sup * 0.5 + x2 * 0.5

            # if the protection leg is under-priced, it means that the hazard rates are under-priced

    # maturity = 2y CDS:
    print("2y CDS with nominal = $1 : premium leg = ",
          premium_leg(spread_2y, N_2y, dntr_2y, intervals_2y, [x1, x1, x2, x2]))
    print("2y CDS with nominal = $1 : protection leg = ",
          protection_leg(r, N_2y, dntr_2y, intervals_2y, [x1, x1, x2, x2]))
    print("hazard rate of 1y ~ 2y = ", x2)
    Prob_1point5y = math.exp(-1 * x2 * 0.5) * Prob_1y
    Prob_2y = math.exp(-1 * x2 * 0.5) * Prob_1point5y
    print("prob survival till 1.5y = ", Prob_1point5y)
    print("prob survival till 2y = ", Prob_2y)

    # Bootstrapping step 3
    print("\n\n# bootstrapping step 3", "\n")
    # assume premiums are paid semi-annually
    # for the maturity = 3 contract
    N_3y = 6
    intervals_3y = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    dntr_3y = [1.03, 1.03, 1.03, 1.03, 1.03, 1.03]
    d = 1
    r = 0.4  # recovery rate
    x_inf = 0.0  # hazard rate trials starting points
    x_sup = 1.0
    x3 = 0.25
    spread_3y = 0.012
    d = 1

    # bisection method to compute hazard rate of stage 3
    while abs(d) > 0.00000001:

        hr = [x1, x1, x2, x2, x3, x3]
        d = protection_leg(r, N_3y, dntr_3y, intervals_3y, hr) - \
            premium_leg(spread_3y, N_3y, dntr_3y, intervals_3y, hr)

        if abs(d) < 0.00000001:
            break

        if d > 0:
            x_sup = x3
            x3 = x_inf * 0.5 + x3 * 0.5
            # if the protection leg is over-priced, it means that the hazard rates are over-priced
        else:
            x_inf = x3
            x3 = x_sup * 0.5 + x3 * 0.5

            # if the protection leg is under-priced, it means that the hazard rates are under-priced

    # maturity = 3y CDS:
    print("3y CDS with nominal = $1 : premium leg = ",
          premium_leg(spread_3y, N_3y, dntr_3y, intervals_3y, [x1, x1, x2, x2, x3, x3]))
    print("3y CDS with nominal = $1 : protection leg = ",
          protection_leg(r, N_3y, dntr_3y, intervals_3y, [x1, x1, x2, x2, x3, x3]))
    print("hazard rate of 2y ~ 3y = ", x3)
    Prob_2point5y = math.exp(-1 * x3 * 0.5) * Prob_2y
    Prob_3y = math.exp(-1 * x3 * 0.5) * Prob_2point5y
    print("prob survival till 2.5y = ", Prob_2point5y)
    print("prob survival till 3y = ", Prob_3y)

    # Bootstrapping step 4
    print("\n\n# bootstrapping step 4", "\n")
    # assume premiums are paid semi-annually
    # for the maturity = 5 contract
    N_5y = 10
    intervals_5y = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    dntr_5y = [1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03]
    d = 1
    r = 0.4  # recovery rate
    x_inf = 0.0  # hazard rate trials starting points
    x_sup = 3.0
    x5 = 0.25
    spread_5y = 0.014

    # bisection method to compute hazard rate of stage 3
    while abs(d) > 0.0000000001:

        hr = [x1, x1, x2, x2, x3, x3, x5, x5, x5, x5]
        d = protection_leg(r, N_5y, dntr_5y, intervals_5y, hr) - \
            premium_leg(spread_5y, N_5y, dntr_5y, intervals_5y, hr)

        if abs(d) < 0.0000000001:
            break

        if d > 0:
            x_sup = x5
            x5 = x_inf * 0.5 + x5 * 0.5
            # if the protection leg is over-priced, it means that the hazard rates are over-priced
        else:
            x_inf = x5
            x5 = x_sup * 0.5 + x5 * 0.5

            # if the protection leg is under-priced, it means that the hazard rates are under-priced

    # maturity = 3y CDS:
    print("5y CDS with nominal = $1 : premium leg = ",
          premium_leg(spread_5y, N_5y, dntr_5y, intervals_5y, [x1, x1, x2, x2, x3, x3, x5, x5, x5, x5]))
    print("5y CDS with nominal = $1 : protection leg = ",
          protection_leg(r, N_5y, dntr_5y, intervals_5y, [x1, x1, x2, x2, x3, x3, x5, x5, x5, x5]))
    print("hazard rate of 3y ~ 5y = ", x5)

    Prob_3point5y = math.exp(-1 * x5 * 0.5) * Prob_3y
    Prob_4y = math.exp(-1 * x5 * 0.5) * Prob_3point5y
    Prob_4point5y = math.exp(-1 * x5 * 0.5) * Prob_4y
    Prob_5y = math.exp(-1 * x5 * 0.5) * Prob_4point5y

    print("prob survival till 3.5y = ", Prob_3point5y)
    print("prob survival till 4y = ", Prob_4y)
    print("prob survival till 4.5y = ", Prob_4point5y)
    print("prob survival till 5y = ", Prob_5y)

    # Visualize the survival curve
    times0 = ['0~0.5', '0.5~1', '1~1.5', '1.5~2', '2~2.5', '2.5~3', '3~3.5', '3.5~4', \
              '4~4.5', '4.5~5']
    times1 = ['0y', '0.5y', '1y', '1.5y', '2y', '2.5y', '3y', '3.5y', '4y', '4.5y', '5y']
    survival_probs = [1, Prob_0point5y, Prob_1y, Prob_1point5y, Prob_2y, Prob_2point5y, Prob_3y, Prob_3point5y, \
                      Prob_4y, Prob_4point5y, Prob_5y]
    hazard_rates = [x1, x1, x2, x2, x3, x3, x5, x5, x5, x5]

    print("answer to question a:")
    print("survival probabilities =", survival_probs)
    print("hazard rates = ", hazard_rates)
    print("the answer to question a is also displayed in two graphs")

    plt.figure()
    plt.subplot(111)
    plt.title('Survival Curve')
    plt.plot(times1, survival_probs, color='tab:blue', marker='o')
    plt.plot(times1, survival_probs, color='black')
    plt.ylabel('survival probability')
    plt.xlabel('time')
    plt.show()

    fig, ax = plt.subplots()
    bar_labels = ['red'] * 10
    ax.bar(times0, hazard_rates, color=bar_labels)

    ax.set_ylabel('hazard rate')
    ax.set_title('Hazard Rates')
    ax.set_xlabel('time intervals : year')
    plt.show()

    # this solves question a
    # we proceed to solve question b
    # we adopt a similar method to compute the fair spread for a 4y CDS
    # Bootstrapping step 4
    print("\n\n# bootstrapping step 4", "\n")
    # assume premiums are paid semi-annually
    # for the maturity = 4 contract
    N_4y = 8
    intervals_4y = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    dntr_4y = [1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03]
    d4 = 1
    r = 0.4  # recovery rate
    s_inf = 0.008  # spread trials starting points
    s_sup = 0.02
    s4 = 0.01
    # this is the spread_4y where we are getting via bisection
    hr = [x1, x1, x2, x2, x3, x3, x5, x5]
    # bisection method to compute hazard rate of stage 3
    while abs(d4) > 0.0000000001:

        d4 = protection_leg(r, N_4y, dntr_4y, intervals_4y, hr) - \
             premium_leg(s4, N_4y, dntr_4y, intervals_4y, hr)

        if abs(d4) < 0.0000000001:
            break

        if d4 < 0:
            s_sup = s4
            s4 = s_inf * 0.5 + s4 * 0.5
            # if the protection leg is under-priced, it means that the spread is over-priced
        else:
            s_inf = s4
            s4 = s_sup * 0.5 + s4 * 0.5

            # if the protection leg is over-priced, it means that the spread is over-priced
    print("4y CDS with nominal = $1 : premium leg = ",
          premium_leg(s4, N_4y, dntr_4y, intervals_4y, hr))
    print("4y CDS with nominal = $1 : protection leg = ",
          protection_leg(r, N_4y, dntr_4y, intervals_4y, hr))
    print("answer to question b:")
    print('the fair price for the 4y CDS spread = ', s4)

    # this solves question b
    # we proceed to solve question c
    # we are paying the original protection seller, say, s5(-1) per year, "-1" means the spread was quoted a year ago.
    # we are receiving the new premium by the amount of s4(0) per year as of today.
    # we should demand [s4(0) - s5(-1)] * RPV01(0)

    s5minus1 = 0.008
    RPV01 = premium_leg(s4, N_4y, dntr_4y, intervals_4y, hr) * pow(s4, -1)
    demand = (s4 - s5minus1) * RPV01
    print("answer to question c:")
    print("The amount I would charge you for transfer of the CDS contract = ", demand)

    # this solves question c
    # we proceed to solve question d
    epsl = 0.000000001

    original_CDS = protection_leg(r, N_4y, dntr_4y, intervals_4y, hr) - \
                   premium_leg(s4 + epsl, N_4y, dntr_4y, intervals_4y, hr) * pow(s4, -1)

    new_CDS = protection_leg(r + epsl, N_4y, dntr_4y, intervals_4y, hr) - \
              premium_leg(s4, N_4y, dntr_4y, intervals_4y, hr) * pow(s4, -1)

    DV01_CDS = (new_CDS - original_CDS) / epsl * -1 * 0.0001
    print("\nanswer to question d")
    print("the DV01 of the 4-y maturity CDS price with repsect to CDS curve = ", DV01_CDS)

    # the following solves question e
    # we suppose there's an upward parallel shift in the interest curve
    # in the following part we compute the CDs price's DV01 with respect to the interest rate curve
    # for the nominal CDS price of $1
    print("\nanswer to question e:")
    dntr_4yn = [1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03, 1.03]

    for i in range(len(dntr_4yn)):
        dntr_4yn[i] = dntr_4yn[i] + epsl

    original_4y = protection_leg(r, N_4y, dntr_4y, intervals_4y, hr) - \
                  premium_leg(s4, N_4y, dntr_4y, intervals_4y, hr) * pow(s4, -1)

    new_4y = protection_leg(r, N_4y, dntr_4yn, intervals_4y, hr) - \
             premium_leg(s4, N_4y, dntr_4yn, intervals_4y, hr) * pow(s4, -1)

    DV01_4y = (new_4y - original_4y) / epsl * -1 * 0.0001

    print("the 4y-maturity CDS price DV01 w.r.t. the interest rate curve =", DV01_4y)

    # question f
    # compute the sensitivity of the CDS price with respect to R
    print("\nanswer to question f:")
    original_price = protection_leg(r, N_4y, dntr_4y, intervals_4y, hr) - \
                     premium_leg(s4, N_4y, dntr_4y, intervals_4y, hr) * pow(s4, -1)

    new_price = protection_leg(r + epsl, N_4y, dntr_4yn, intervals_4y, hr) - \
                premium_leg(s4, N_4y, dntr_4yn, intervals_4y, hr) * pow(s4, -1)

    sensitivity_r = (new_price - original_price) / epsl
    print("the sensitivity of the 4 year CDS price with repsect to R = ", sensitivity_r)
