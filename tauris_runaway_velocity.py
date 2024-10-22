G = 6.67430e-11  # m^3 kg^-1 s^-2

# Define the runaway velocity function (as provided)
def runaway_velocity(M1, M_CO, M2, r, omega, theta, phi, v_im):

    # Calculate m_shell (mass of the ejected shell)
    m_shell = (M1 - M_CO) / M_CO
    m_2 = M2 / M_CO
    
    v = np.sqrt(G * (M1 + M2) / r)
    
    m_tilde = (1 + m_2) / (1 + m_shell + m_2)
    
    P = 1 - 2 * m_tilde + (omega ** 2 / v ** 2) + (v_im ** 2 / v ** 2) + (2 * omega / v ** 2) * (v * np.cos(theta) - v_im * np.sin(theta) * np.cos(phi))
    
    Q = 1 + (P / m_tilde) - ((omega * np.sin(theta) * np.cos(phi) - v_im) ** 2) / (m_tilde * v ** 2)
    
    R = ((np.sqrt(np.abs(P)) / (m_tilde * v)) * (omega * np.sin(theta) * np.cos(phi) - v_im) - (P / m_tilde) - 1) * (1 + m_2) / m_2
    
    S = (1 + (P * (Q + 1) / m_tilde)) * (1 + m_2) / m_2
    
    v_x = -(omega * np.cos(theta)) / (m_2 * R) - ((1 / (m_2 * R)) + (1 + m_shell) / (1 + m_shell + m_2)) * v
    v_y = (omega * np.sin(theta) * np.cos(phi)) / (m_2 * S) + (1 - 1 / (m_2 * S)) * v_im - (Q * np.sqrt(np.abs(P)) / (m_2 * S)) * v
    v_z = -(omega * np.sin(theta) * np.sin(phi)) / (m_2 * R)
    
    return v_x, v_y, v_z

def monte_carlo_runaway(df, iterations=10000, sigma = 190e3 / np.sqrt(2)):
    velocities = []
    prob = []
    for index, row in df.iterrows():

        M1 = row['m1_prev'] * 1e30  #primary mass before supernova
        M_CO = row['m1'] * 1e30  #ompact object mass after supernova
        M2 = row['m2'] * 1e30  #companion mass
        r = np.abs(row['a']) * 1.496e11  #orbital separation (AU to meters)
        r2 = row['r2'] * 6.96e8  #companion radius
        probability = row['prob']

        A = 1.0664e3 - 2.6154e2 * M2 / M_CO
        eta1 = -1.85690 + 7.69070e-3 * M2 / M_CO
        v_im = A * np.power(r / r2, eta1)

        for _ in range(iterations):
            #omega = np.random.normal(vkick, sigma)  # Using mu=190 km/s and sigma from user input
            omega = stats.maxwell.rvs(scale=sigma)
            theta = np.random.uniform(0, np.pi * 2)  # theta in range [0, 2pi]
            phi = np.random.uniform(0, np.pi)  # phi in range [0, pi]

            v_x, v_y, v_z = runaway_velocity(M1, M_CO, M2, r, omega, theta, phi, v_im)
            total_velocity = np.sqrt(v_x**2 + v_y**2 + v_z**2)/1e3
            
            velocities.append(total_velocity)
            prob.append(probability)
    
    return velocities, prob
