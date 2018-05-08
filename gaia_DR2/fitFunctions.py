def gauss(x, *p):
   # A, mu, sigma = p
   # return A*np.exp(-(x-mu)**2/(2.*sigma**2))
    A, mu, sigma = p
    return A*np.exp(-(x-0)**2/(2.*sigma**2))

def fit_f_v_z(hist, bin_edges, name):
    # Define some test data which is close to Gaussian
    
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    
    # p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
    p0 = [1., 0., 1.]
    
    coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)
    
    # Get the fitted curve
    hist_fit = gauss(bin_centres, *coeff)
    
    plt.plot(bin_centres, hist, label='Test data')
    plt.plot(bin_centres, hist_fit, label='Fitted data')
    plt.xlabel("velocity [km/s]")
    plt.ylabel("count")
    plt.title(name)
    
    # Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
    print 'Fitted mean = ', coeff[1]
    print 'Fitted height = ', coeff[0]
    print 'Fitted standard deviation = ', coeff[2]
    #plt.show()
    plt.savefig("output/f_v_z_{}.png".format(name))

    return coeff