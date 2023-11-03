#include "../binary_c.h"
No_empty_translation_unit_warning;

#include "supernovae.h"

void tauris_takens_orbit(struct stardata_t * const stardata,
                         struct stardata_t * const pre_explosion_stardata,
                         struct star_t * const pre_explosion_star,
                         struct star_t * const exploder,
                         struct star_t * const companion,
                         struct kick_system_t * const kick_system,
                         const int kick_distribution MAYBE_UNUSED,
                         const double kick_dispersion MAYBE_UNUSED,
                         const double kick_companion MAYBE_UNUSED)
{
    /*
     * Overwrite the previously evaluated kick_system->new_orbital_speed_squared
     *
     * if entering this branch, calculate the final velocity
     * according to Takens & Tauris (1998). This assumes the orbit
     * pre-SN explosion is circular but accounts for
     * the momentum transfer from the ejecta shell to the secondary,
     * giving more realistic final velocities.
     *
     * In this section the same notation as in
     * Tauris & Takens (1998) is used, except:
     *
     * w  ->  kick_system->kick_speed
     * m2 ->  m2i to avoid overlaps with previous notation
     *
     * N.B: m2i is the initial mass in units of the remnant mass
     */

    /*
     * small m are masses in units of the remnant mass
     */
    double v_im=0.0,mtilde,m2f,msh,m2i,Mremnant,r2;

    /* parameter initialization */
    Mremnant = Max(1e-50,exploder->mass);
    m2i = Max(1e-50,companion->mass / Mremnant);
    r2 = companion->radius;
    Dprint("MremnaNt %g, m2i = %g, r2 = %g\n",Mremnant,m2i,r2);

    /*
     * Fraction of mass stripped from the secondary by the SN ejecta.
     * Note that we do not distinguish between methods of stripping:
     * it could be shock heated ablation or direct stripping.
     */
    double Fstrip = 0.0;

    /*
     * Fraction of the ejecta that is accreted by the secondary
     */
    double Facc = 0.0;

    /*
     * Efficiency of impact //0.5;// value from TT98
     */
    double eta=1.0;
    double xcrit=1.0;

    /*
     * typical SN ejecta velocity in km/s
     */
    double v_ejecta = supernova_ejecta_velocity(exploder);

    /*
     * Do we want to kick or strip the companion?
     */
    if(Kick_companion_with(SN_IMPULSE_LIU2015))
    {
        if(ON_EITHER_MAIN_SEQUENCE(companion->stellar_type))
        {
            /*
             * Use Liu et al. (2015, A&A 584, A11) to determine
             * how much material is lost from a main-sequence star.
             *
             * The following expressions interpolate their
             * fitting coefficients from their 0.9 and 3.5Msun
             * main-sequence secondary models to assess the total mass
             * stripped and ablated from the SN ejecta hitting the secondary
             * (they show deltaM2/M2 which is what we want).
             *
             * See Fig. 6, leftmost panel.
             */
            double x = stardata->common.orbit.separation / companion->radius;
            double A = Max(0.0, 5.69860 - 1.58730 * companion->mass);
            double eta1 = -3.4442 + 2.26920e-01 * companion->mass;
            Fstrip = A * pow(x,eta1);

            /*
             * limit to their ~ max of 10%
             */
            Clamp(Fstrip, 0.0, 0.1);

            /*
             * Below is a similar fit for accretion onto the companion
             * star (the right panel of Fig. 6) but note that their
             * fitting function gives the amount of accreted material
             * in Msun. We want the fraction of the companion stripped,
             * i.e. Delta M / companion->mass,
             * as Liu+ plot in the leftmost panel of Fig. 6.
             */
            A = (4.03460e-03)+(1.96150e-03)*companion->mass;
            eta1 = (-7.81540e-01)+(-1.13850e-01)*companion->mass;
            Facc = A * pow(x,eta1);

            /*
             * Limit to 1e-2Msun
             */
            Clamp(Facc, 0.0, 1e-2);

            /* normalize to M2 */
            Facc /= companion->mass;
        }
        else if(GIANT_LIKE_STAR(companion->stellar_type))
        {
            /*
             * Giants from Hirai+ (2014, ApJ 792,66)
             * Eq.4 and Figure 16.
             * This formula gives Mej/M2 (i.e. it is
             * already normalized to m2i) as a function of
             * a/amin where amin is the minimum radius for
             * RLOF. This factor is 1.0/rl(q) where q=M2/M1
             * and star 1 is the exploder.
             *
             * TODO use Pan et al 2012 (ApJ 750,151) and
             * Pan et al 2013 (ApJ 773,49) for giants.
             */
            double A = 0.26;
            double eta1 = -4.3;
            double q = companion->mass / pre_explosion_star->mass;
            double x = 1.0 / rl(q);
            Fstrip = A * pow(x,eta1);

            /*
             * At a/amin = 1 they find a maximum strip fraction
             * of ~ 0.25.
             */
            Clamp(Fstrip, 0.0, 0.25);
        }
    }

    Dprint("Fstrip = %g : Facc = %g : companion type %d\n",
           Fstrip,
           Facc,
           companion->stellar_type);

    // mass of the ejected shell in units of the remnant mass
    msh = pre_explosion_star->mass / Mremnant - 1.0 ;

    m2f = m2i * (1.0 - Fstrip + Facc); // m2f == m2i if Fstrip + Facc = 0

    mtilde = (1.0 + m2f) / (1.0 + msh + m2i);
    Dprint("mtilde = %g from m2f = %g, msh = %g, m2i = %g\n",mtilde,m2f,msh,m2i);

    /*
     * The orbital velocity prior to explosion
     */
    double orbv = stardata->model.sgl == TRUE ? 1e-50 :
        Max(1e-50,orbital_velocity(pre_explosion_stardata));
    Dprint("orbv = %g\n",orbv);
    SNprint("SN orbv2 from %g %g %g\n",
            pre_explosion_stardata->star[0].mass,
            pre_explosion_stardata->star[1].mass,
            pre_explosion_stardata->common.orbit.separation);

    SNprint("SN orbv = %g\n",orbv);

    if(Kick_companion_with(SN_IMPULSE_LIU2015))
    {
        if(SN_TYPE_IS_IA(exploder->SN_type))
        {
            /*
             * Implement Liu+ (2012 A&A 548, A2)
             */
            Dprint("Kick companion with SN_IMPULSE_LIU2015 (2012 for SNeIa)");
            v_im = 0.0; // Eq. 3
        }
        else if(SN_TYPE_IS_CORE_COLLAPSE(exploder->SN_type))
        {
            /*
             * Liu+2015 give a fitting formula for 0.9 and 3.5Msun
             * main-sequence companions with a core collapse supernovae.
             *
             * Interpolate the coefficients as a function of mass.
             *
             * NB Limit the mass (set in Mcomp) to the range fitted
             * by Liu+2015, otherwise the expression can become negative,
             * which is clearly incorrect.
             */
            const double Mcomp = Limit_range(companion->mass,0.9,3.5);
            Dprint("Kick companion with SN_IMPULSE_LIU2015 (2015 for CCSNe): companion->mass = %g : use Mcomp = %g",
                   companion->mass,
                   Mcomp);
            const double A = 1.0664e3 - 2.6154e2 * Mcomp;
            const double eta1 =  -1.85690  + 7.69070e-3 * Mcomp;
            v_im = A * pow(stardata->common.orbit.separation/r2,eta1);
        }
        else
        {
            v_im = 0.0;
        }
        Dprint("Impact speed %g\n",v_im);
    }
    else if(Kick_companion_with(SN_IMPULSE_WHEELER1975) &&
            exploder->SN_type != SN_WDKICK)
    {
        /*
         * use fitting formula from Wheeler, Lecar & McKee (1975, Eq. 26)
         * Eq. 5 in Tauris & Takens (1998, A&A 330, 1047).
         */
        Dprint("Kick companion with SN_IMPULSE_WHEELER1975");
        const double v_esc = sqrt(GMRKM*m2i*Mremnant/(xcrit*r2));
        const double rfac = 0.5*r2/stardata->common.orbit.separation;

        v_im = eta * v_ejecta * Pow2(rfac) * (msh/m2i) * Pow2(xcrit) *
            (1.0+log(2.0*v_ejecta/v_esc))/(1.0 - Fstrip + Facc);

        if(!Fequal(xcrit,1.0))
        {
            SNprint("MATHIEU_WARNING: xcrit=%e!=1 => you need to change the formula for the impact velocity!\n",
                    xcrit);
        }
    }
    else
    {
        v_im = 0.0;
    }

    /*
     * The estimates of v_im are very uncertain in both cases. Let's force it to be
     * greater than the escape velocity from the non-exploding star, since if the ejecta
     * were free-falling we would have this velocity (so it's a good lower boundary)
     */
    if(0 && v_im > TINY)
    {
        v_im = Max(v_im, sqrt(GMRKM*m2i*Mremnant/r2));
    }
    Dprint("VIM %g\n",v_im);

    /*
     * P,Q,R,S are defined to simplify notation
     */
    double P,Q,R,S;

    /*
     * P -- Eq.44 in Tauris & Takens 98
     *
     * NB if P<0 then the system remains bound,
     * and v_cm is given by TT98's Eq.6 (see
     * their Sec. 3.3 on p1058). At present, we
     * then ignore the kick on the companion
     * and simply return (to use the hn_squared of Hurley+2002).
     *
     * This is incorrect and requires a fix.
     */
    P = 1.0
        -
        2.0 * mtilde
        +
        Pow2(kick_system->kick_speed) / Pow2(orbv)
        +
        Pow2(v_im)/Pow2(orbv)
        +
        2.0 * kick_system->kick_speed / Pow2(orbv)*
        (orbv*kick_system->cosomega - v_im*kick_system->sinomega*kick_system->cosphi);

    Dprint("P = 1.0 - %g + %g + %g(impact) + %g * %g = %g\n",
           2.0 * mtilde,
           Pow2(kick_system->kick_speed) / Pow2(orbv),
           Pow2(v_im)/Pow2(orbv),
           2.0 * kick_system->kick_speed / Pow2(orbv),
           (orbv*kick_system->cosomega - v_im*kick_system->sinomega*kick_system->cosphi),
           P);

    double v2x,v2y,v2z;
    double vrem_x,vrem_y,vrem_z;
    double vnx,vny,vnz;

    if(P<0.0)
    {
        Dprint("P<0 : system is bound\n");
        return;
    }
    else
    {
        /*
         * Q -- Eq.45 in Tauris & Takens 98
         */
        Q = 1.0 + P / mtilde
            -
            Pow2(kick_system->kick_speed * kick_system->sinomega * kick_system->cosphi - v_im)
            /
            (mtilde*Pow2(orbv));
        Dprint("Q = 1 + %g / %g - (%g * %g * %g - %g)^2 / (%g * %g^2) = %g\n",
               P, mtilde,
               kick_system->kick_speed,

               kick_system->sinomega,
               kick_system->cosphi,
               v_im,
               mtilde,
               orbv,
               Q);

        /*
         * R -- Eq.46 in Tauris & Takens 98
         */
        double Rsub =  (kick_system->kick_speed * kick_system->sinomega * kick_system->cosphi - v_im);

        Dprint("Rsub = %g\n",Rsub);

        R = (
            Eval_if_nonzero(Rsub, sqrt(P) / (mtilde*orbv) * Rsub)
            -
            P / mtilde
            -
            1.0
            )
            *
            (
                1.0/m2f
                +
                1.0
                );

        Dprint("R = %g from P=%g m2f=%g mtilde=%g orbv=%g kick_speed=%g\n",
               R,
               P,
               m2f,
               mtilde,
               orbv,
               kick_system->kick_speed);
        Dprint(" sinomega=%g cosphi=%g v_im=%g m2f=%g",
               kick_system->sinomega,
               kick_system->cosphi,
               v_im,
               m2f);

        /*
         * S -- Eq.47 in Tauris & Takens 98
         */
        S = (
            1.0
            +
            P * (Q + 1.0) / mtilde
            )
            *
            (
                1.0 / m2f
                +
                1.0
                );
        Dprint("S = %g from P=%g Q=%g mtilde=%g m2f=%g\n",
               S,P,Q,mtilde,m2f);

        /*
         * The following are at infinity and in the reference
         * frame of the total center of mass
         * (including ejecta masses), i.e. the original frame
         */

        /*
         * Remnant velocity components (Eqs. 51-53 in Tauris & Takens 1998)
         */

        /*
         * companion velocity (Eqs. 54-56 in Tauris & Takens 1998)
         */
        vrem_x =
            kick_system->kick_speed*
            kick_system->cosomega*
            (1.0 + 1.0 / R)
            +
            orbv *
            (1.0/R + m2i/(1.0+msh+m2i));
        SNprint("SN vrem_x = %g * %g * %g + %g * (%g + %g) = %g\n",
                kick_system->kick_speed,
                kick_system->cosomega,
                (1.0 + 1.0 / R),
                orbv,
                1.0/R,
                m2i/(1.0+msh+m2i),
                vrem_x
            );
        Dprint("VREMX = %g * %g * %g + %g * (%g + %g) = %g\n",
               kick_system->kick_speed,
               kick_system->cosomega,
               (1.0 + 1.0 / R),
               +orbv,
               Eval_if_nonzero(R,1.0/R),
               + m2i/(1.0+msh+m2i),
               vrem_x
            );

        vrem_y =
            kick_system->kick_speed*
            kick_system->sinomega*
            kick_system->cosphi*
            (1.0 - 1.0 / S)
            +
            v_im / S
            +
            Eval_if_nonzero(Q, orbv*Q*sqrt(P) / S);

        vrem_z =
            kick_system->kick_speed*
            kick_system->sinomega*
            kick_system->sinphi*
            (1.0 + 1.0/R);

        if(m2f > TINY)
        {
            /* Eq. 54 */
            v2x = -kick_system->kick_speed*kick_system->cosomega/(m2f*R)
                -
                orbv*(1.0/(m2f * R) + (1.0 + msh)/(1.0 + msh + m2i));

            SNprint("v2x = - %g * %g / (%g * %g) - %g * (1/(%g * %g) + (1+%g)/(1+%g+%g)) = %g - %g * %g = %g - %g = %g\n",
                    kick_system->kick_speed,
                    kick_system->cosomega,
                    m2f,
                    R,
                    orbv,
                    m2f,
                    R,
                    msh,
                    msh,
                    m2i,
                    -kick_system->kick_speed*kick_system->cosomega/(m2f*R),
                    orbv,
                    (1.0/(m2f * R) + (1.0 + msh)/(1.0 + msh + m2i)),
                    -kick_system->kick_speed*kick_system->cosomega/(m2f*R),
                    orbv*(1.0/(m2f * R) + (1.0 + msh)/(1.0 + msh + m2i)),
                    v2x);

            /* Eq. 55 */
            v2y = kick_system->kick_speed*kick_system->sinomega*kick_system->cosphi/(m2f*S)
                +
                v_im*(1.0 - 1.0/(m2f*S))
                -
                Eval_if_nonzero(Q, Q*orbv*sqrt(P)/(m2f*S));

            SNprint("v2y = %g * %g * %g / (%g * %g) + (1-1/(%g * %g))*%g - %g * %g * %g / (%g * %g) = %g\n",
                    kick_system->kick_speed,
                    kick_system->sinomega,
                    kick_system->cosphi,
                    m2f,
                    S,
                    v_im,
                    m2f,
                    S,
                    orbv,
                    Q,
                    Eval_if_nonzero(Q, sqrt(P)),
                    m2f,
                    S,
                    v2y);

            v2z = -kick_system->kick_speed*kick_system->sinomega*kick_system->sinphi/(m2f*R);
        }
        else
        {
            v2x = v2y = v2z = 0.0;
        }

        /*
         * Calculate the new relative velocity vector
         */
        vnx = vrem_x - v2x;
        SNprint("SN vnx=%g from vrem_x=%g v2x=%g\n",
                vnx,vrem_x,v2x);
        vny = vrem_y - v2y;
        vnz = vrem_z - v2z;
        Dprint ("vnx = %g - %g = %g\n",vrem_x,v2x,vnx);
        Dprint ("vny = %g - %g = %g\n",vrem_y,v2y,vny);
        Dprint ("vnz = %g - %g = %g\n",vrem_z,v2z,vnz);


        /* now update the quantities used to evaluate the new orbit */
        kick_system->new_orbital_speed_squared = Pythag3_squared(vnx,vny,vnz);
        kick_system->new_orbital_speed = sqrt(kick_system->new_orbital_speed_squared);

        SNprint("SN vnx=%g vny=%g vnz=%g : new orbital speed %g\n",
                vnx,vny,vnz,kick_system->new_orbital_speed_squared);

    }

    /*
     * Set the stripped/ablated secondary mass
     *
     * recall: m2f was in units of the remnant mass, but m2
     * is instead the new mass in solar units of the secondary!!
     */
    companion->mass = m2f * Mremnant;
    if(companion->mass < TINY)
    {
        /* completely disrupted */
        stellar_structure_make_massless_remnant(stardata,
                                                companion);
    }

#ifdef NUCSYN
    companion->dm_companion_SN = m2f - m2i;
#endif //NUCSYN

    /*
     * and save the velocity in the variables for runaway velocity
     */
#ifdef RUNAWAY_STARS
    exploder->v_Rwx += vrem_x;
    exploder->v_Rwy += vrem_y;
    exploder->v_Rwz += vrem_z;
    exploder->v_Rw = Pythag3(vrem_x,vrem_y,vrem_z);

    companion->v_Rwx = v2x;
    companion->v_Rwy = v2y;
    companion->v_Rwz = v2z;
    companion->v_Rw = Pythag3(v2x,v2y,v2z);

    Dprint("vremx=%e,vremy=%e,vremz=%e\n",vrem_x,vrem_y,vrem_z);
    Dprint("----------------------------------------------------\n");
    Dprint("vnx=%e,vny=%e,vnz=%e\n",vnx,vny,vnz);
    Dprint("----------------------------------------------------\n");
    Dprint("vn2=%e,m1n=%e,m2f=%e\n",kick_system->new_orbital_speed_squared,Mremnant,m2f);
    Dprint("----------------------------------------------------\n");
    Dprint("v_Rw1=%e, v_Rw_exploding=%e, v_Rw2=%e\n",
           stardata->star[0].v_Rw,
           exploder->v_Rw,
           stardata->star[1].v_Rw);
    Dprint("----------------------------------------------------\n");
    Dprint("v2x=%e,v2y=%e,v2z=%e\n",v2x,v2y,v2z);
    Dprint("----------------------------------------------------\n");
    Dprint("v_Rwx1=%e, v_Rwy1=%e, v_Rwz1=%e\n",
           stardata->star[0].v_Rwx,
           stardata->star[0].v_Rwy,
           stardata->star[0].v_Rwz);
    Dprint("v_Rwx2=%e, v_Rwy2=%e, v_Rwz2=%e\n",
           stardata->star[1].v_Rwx,
           stardata->star[1].v_Rwy,
           stardata->star[1].v_Rwz);
    Dprint("----------------------------------------------------\n");
#endif


    /*
     * r along y-direction by definition of the ref. frame
     * (both in BSE and TT98)
     *
     * hn_squared is the squared modulus of the new
     * specific orbital angular momentum:
     *
     * hn_squared=|hn|^2
     *
     * with: hn = r x vn = (r*vnz, 0, r*vnx) => hn_squared = r^2(vnx^2+vnz^2)
     */
    kick_system->hn_squared = Pow2(kick_system->separation) *
        Pythag2_squared(vnx,vnz);

    Dprint("hn^2 = %g^2 * (vnx^2=%g + vnz^2=%g) = %g (Pythag2_squared = %g)\n",
           kick_system->separation,
           vnx*vnx,
           vnz*vnz,
           kick_system->hn_squared,
           Pythag2_squared(vnx,vnz)
        );

    kick_system->hn = sqrt(kick_system->hn_squared);

}
