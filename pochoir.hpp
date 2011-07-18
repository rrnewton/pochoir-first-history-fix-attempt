/*
 **********************************************************************************
 *  Copyright (C) 2010-2011  Massachusetts Institute of Technology
 *  Copyright (C) 2010-2011  Yuan Tang <yuantang@csail.mit.edu>
 * 		                     Charles E. Leiserson <cel@mit.edu>
 * 	 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Suggestsions:                  yuantang@csail.mit.edu
 *   Bugs:                          yuantang@csail.mit.edu
 *
 *********************************************************************************
 */

#ifndef EXPR_STENCIL_HPP
#define EXPR_STENCIL_HPP

#include "pochoir_common.hpp"
#include "pochoir_types.hpp"
#include "pochoir_kernel.hpp"
#include "pochoir_array.hpp"
/* assuming there won't be more than 10 Pochoir_Array in one Pochoir object! */
#define ARRAY_SIZE 10

template <int N_RANK>
class Pochoir {
    private:
        int slope_[N_RANK];
        grid_info<N_RANK> logic_grid_;
        grid_info<N_RANK> phys_grid_;
        int time_shift_;
        int toggle_;
        int timestep_;
        bool regArrayFlag, regLogicDomainFlag, regPhysDomainFlag, regShapeFlag;
        void checkFlag(bool flag, char const * str);
        void checkFlags(void);
        template <typename T_Array>
        void getPhysDomainFromArray(T_Array & arr);
        template <typename T_Array>
        void cmpPhysDomainFromArray(T_Array & arr);
        void Register_Shape(Pochoir_Shape<N_RANK> * shape, int N_SIZE);
        Pochoir_Shape<N_RANK> * shape_;
        int shape_size_;
        int num_arr_;
        int arr_type_size_;
        int sz_pgk_;
        int lcm_unroll_;
        /* assuming that the number of distinct sub-regions is less than 10 */
        Pochoir_Guard_Kernel<N_RANK> * pgk_;
        Pochoir_Obase_Guard_Kernel<N_RANK> * opgk_;

        /* Private Register Kernel Function */
        template <typename K>
        void reg_kernel(int pt, K k);
        template <typename K, typename ... KS>
        void reg_kernel(int pt, K k, KS ... ks);

    public:
    template <typename ... KS>
    void Register_Kernel(typename Pochoir_Types<N_RANK>::T_Guard g, KS ... ks);
    void Register_Obase_Kernel(typename Pochoir_Types<N_RANK>::T_Guard g, int unroll, typename Pochoir_Types<N_RANK>::T_Obase_Kernel k, typename Pochoir_Types<N_RANK>::T_Obase_Kernel cond_k, typename Pochoir_Types<N_RANK>::T_Obase_Kernel bk, typename Pochoir_Types<N_RANK>::T_Obase_Kernel cond_bk);
    // get slope(s)
    int slope(int const _idx) { return slope_[_idx]; }
    Pochoir() {
        for (int i = 0; i < N_RANK; ++i) {
            slope_[i] = 0;
            logic_grid_.x0[i] = logic_grid_.x1[i] = logic_grid_.dx0[i] = logic_grid_.dx1[i] = 0;
            phys_grid_.x0[i] = phys_grid_.x1[i] = phys_grid_.dx0[i] = phys_grid_.dx1[i] = 0;
        }
        shape_ = NULL; shape_size_ = 0; time_shift_ = 0; toggle_ = 0;
        timestep_ = 0;
        regArrayFlag = regLogicDomainFlag = regPhysDomainFlag = regShapeFlag = false;
        num_arr_ = 0;
        arr_type_size_ = 0;
        sz_pgk_ = 0;
        pgk_ = NULL; opgk_ = NULL; 
        lcm_unroll_ = 1;
    }

    /* currently, we just compute the slope[] out of the shape[] */
    /* We get the grid_info out of arrayInUse */
    template <typename A>
    void Register_Array(A & a);
    template <typename A, typename ... AS>
    void Register_Array(A & a, AS ... as);

    /* We should still keep the Register_Domain for zero-padding!!! */
    template <typename D>
    void Register_Domain(D const & d);
    template <typename D, typename ... DS>
    void Register_Domain(D const & d, DS ... ds);
    /* register boundary value function with corresponding Pochoir_Array object directly */
    template <typename T_Array, typename RET>
    void registerBoundaryFn(T_Array & arr, RET (*_bv)(T_Array &, int, int, int)) {
        arr.Register_Boundary(_bv);
        Register_Array(arr);
    } 
    grid_info<N_RANK> get_phys_grid(void);

    /* Executable Spec */
    void Run(int timestep);
    void Run_Obase(int timestep);
    /* obase for zero-padded region */
    template <typename F>
    void Run_Obase(int timestep, F const & f);
    /* obase for interior and ExecSpec for boundary */
    template <typename F, typename BF>
    void Run_Obase(int timestep, F const & f, BF const & bf);
};

template <int N_RANK> template <typename K>
void Pochoir<N_RANK>::reg_kernel(int pt, K k) {
    pgk_[sz_pgk_].kernel_[pt] = k; 
    Register_Shape(k.Get_Shape(), k.Get_Shape_Size());
}

template <int N_RANK> template <typename K, typename ... KS>
void Pochoir<N_RANK>::reg_kernel(int pt, K k, KS ... ks) {
    pgk_[sz_pgk_].kernel_[pt] = k;
    Register_Shape(k.Get_Shape(), k.Get_Shape_Size());
    reg_kernel(pt+1, ks ...);
}

template <int N_RANK> template <typename ... KS>
void Pochoir<N_RANK>::Register_Kernel(typename Pochoir_Types<N_RANK>::T_Guard g, KS ... ks) {
    int l_size = sizeof...(KS);
//    typedef typename Pochoir_Types<N_RANK>::T_Kernel T_Kernel;
    typedef Pochoir_Kernel<N_RANK> T_Kernel;
    if (pgk_ == NULL) {
        pgk_ = new Pochoir_Guard_Kernel<N_RANK>[ARRAY_SIZE];
        sz_pgk_ = 0;
        assert(lcm_unroll_ == 1);
    }
    assert(sz_pgk_ < ARRAY_SIZE);
    if (sz_pgk_ >= ARRAY_SIZE) {
        printf("Pochoir Error: Register_Kernel > %d\n", sz_pgk_);
        exit(1);
    }
    pgk_[sz_pgk_].unroll_ = l_size;
    pgk_[sz_pgk_].pointer_ = 0;
    pgk_[sz_pgk_].guard_ = g;
    pgk_[sz_pgk_].kernel_ = (T_Kernel *) calloc(l_size, sizeof(T_Kernel));
    reg_kernel(0, ks ...);
    lcm_unroll_ = lcm(lcm_unroll_, l_size);
    ++sz_pgk_;
}

template <int N_RANK> 
void Pochoir<N_RANK>::Register_Obase_Kernel(typename Pochoir_Types<N_RANK>::T_Guard g, int unroll, typename Pochoir_Types<N_RANK>::T_Obase_Kernel k, typename Pochoir_Types<N_RANK>::T_Obase_Kernel cond_k, typename Pochoir_Types<N_RANK>::T_Obase_Kernel bk, typename Pochoir_Types<N_RANK>::T_Obase_Kernel cond_bk) {
    typedef typename Pochoir_Types<N_RANK>::T_Obase_Kernel T_Kernel;
    if (opgk_ == NULL) {
        opgk_ = new Pochoir_Obase_Guard_Kernel<N_RANK>[ARRAY_SIZE];
        sz_pgk_ = 0;
        assert(lcm_unroll_ == 1);
    }
    assert(sz_pgk_ < ARRAY_SIZE);
    if (sz_pgk_ >= ARRAY_SIZE) {
        printf("Pochoir Error: Register_Kernel > %d\n", sz_pgk_);
        exit(1);
    }
    opgk_[sz_pgk_].guard_ = g;
    opgk_[sz_pgk_].unroll_ = unroll;
    opgk_[sz_pgk_].kernel_ = k;
    opgk_[sz_pgk_].cond_kernel_ = cond_k;
    opgk_[sz_pgk_].bkernel_ = bk;
    opgk_[sz_pgk_].cond_bkernel_ = cond_bk;
    lcm_unroll_ = lcm(lcm_unroll_, unroll);
    ++sz_pgk_;
}

template <int N_RANK>
void Pochoir<N_RANK>::checkFlag(bool flag, char const * str) {
    if (!flag) {
        printf("\nPochoir registration error:\n");
        printf("You forgot to register %s.\n", str);
        exit(1);
    }
}

template <int N_RANK>
void Pochoir<N_RANK>::checkFlags(void) {
    checkFlag(regArrayFlag, "Pochoir array");
    checkFlag(regLogicDomainFlag, "Logic Domain");
    checkFlag(regPhysDomainFlag, "Physical Domain");
    checkFlag(regShapeFlag, "Shape");
    return;
}

template <int N_RANK> template <typename T_Array> 
void Pochoir<N_RANK>::getPhysDomainFromArray(T_Array & arr) {
    /* get the physical grid */
    for (int i = 0; i < N_RANK; ++i) {
        phys_grid_.x0[i] = 0; phys_grid_.x1[i] = arr.size(i);
        /* if logic domain is not set, let's set it the same as physical grid */
        if (!regLogicDomainFlag) {
            logic_grid_.x0[i] = 0; logic_grid_.x1[i] = arr.size(i);
        }
    }

    regPhysDomainFlag = true;
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename T_Array> 
void Pochoir<N_RANK>::cmpPhysDomainFromArray(T_Array & arr) {
    /* check the consistency of all engaged Pochoir_Array */
    for (int j = 0; j < N_RANK; ++j) {
        if (arr.size(j) != phys_grid_.x1[j]) {
            printf("Pochoir array size mismatch error:\n");
            printf("Registered Pochoir arrays have different sizes!\n");
            exit(1);
        }
    }
}

template <int N_RANK> template <typename A>
void Pochoir<N_RANK>::Register_Array(A & a) {
    if (!regShapeFlag) {
        printf("Please register Shape before register Array!\n");
        exit(1);
    }

    if (num_arr_ == 0) {
        arr_type_size_ = sizeof(A);
        ++num_arr_;
    } 
    if (!regPhysDomainFlag) {
        getPhysDomainFromArray(a);
    } else {
        cmpPhysDomainFromArray(a);
    }
    a.Register_Shape(shape_, shape_size_);
#if 0
    arr.set_slope(slope_);
    arr.set_toggle(toggle_);
    arr.alloc_mem();
#endif
    regArrayFlag = true;
}

template <int N_RANK> template <typename A, typename ... AS>
void Pochoir<N_RANK>::Register_Array(A & a, AS ... as) {
    if (!regShapeFlag) {
        printf("Please register Shape before register Array!\n");
        exit(1);
    }

    if (num_arr_ == 0) {
        arr_type_size_ = sizeof(A);
        ++num_arr_;
    } 
    if (!regPhysDomainFlag) {
        getPhysDomainFromArray(a);
    } else {
        cmpPhysDomainFromArray(a);
    }
    a.Register_Shape(shape_, shape_size_);
#if 0
    arr.set_slope(slope_);
    arr.set_toggle(toggle_);
    arr.alloc_mem();
#endif
    Register_Array(as ...);
    regArrayFlag = true;
}

template <int N_RANK> 
void Pochoir<N_RANK>::Register_Shape(Pochoir_Shape<N_RANK> * shape, int N_SIZE) {
    Pochoir_Shape<N_RANK> * l_shape;
    if (!regShapeFlag) {
        regShapeFlag = true;
        l_shape = new Pochoir_Shape<N_RANK>[N_SIZE];
        shape_ = l_shape;
    } else {
        l_shape = new Pochoir_Shape<N_RANK>[N_SIZE + shape_size_] ;
        /* copy in the old shape value */
        for (int i = 0; i < shape_size_; ++i)  {
            for (int r = 0; r < N_RANK+1; ++r)  {
                l_shape[i].shift[r]  = shape_[i].shift[r];
            }
        }
        delete (shape_);
        shape_ = l_shape;
    }
    /* currently we just get the slope_[]  and toggle_ out of the shape[]  */
    int l_min_time_shift=0, l_max_time_shift=0, depth=0;
    for (int i = 0; i < N_SIZE; ++i)  {
        if (shape[i] .shift[0] < l_min_time_shift)
            l_min_time_shift = shape[i].shift[0];
        if (shape[i].shift[0] > l_max_time_shift)
            l_max_time_shift = shape[i].shift[0];
    }
    depth = l_max_time_shift - l_min_time_shift;
    time_shift_ = max(time_shift_, 0 - l_min_time_shift);
    toggle_ = max(toggle_, depth + 1);
    for (int i = 0; i < N_SIZE; ++i) {
        for (int r = 0; r < N_RANK+1; ++r) {
            slope_[N_RANK-r] = (r > 0) ? max(slope_[N_RANK-r], abs((int)ceil((float)shape_[i].shift[r]/(l_max_time_shift - shape_[i].shift[0])))) : 0;
            shape_[shape_size_ + i].shift[r] = (r > 0) ? shape[i].shift[r] : shape[i].shift[r] + time_shift_;
        }
    }
    shape_size_ += N_SIZE;
#if DEBUG 
    printf("time_shift_ = %d, toggle = %d\n", time_shift_, toggle_);
    for (int r = 0; r < N_RANK; ++r) {
        printf("slope[%d] = %d, ", r, slope_[r]);
    }
    printf("\n");
#endif
}

template <int N_RANK> template <typename D>
void Pochoir<N_RANK>::Register_Domain(D const & d) {
    logic_grid_.x0[0] = d.first();
    logic_grid_.x1[0] = d.first() + d.size();
    regLogicDomainFlag = true;
}

template <int N_RANK> template <typename D, typename ... DS>
void Pochoir<N_RANK>::Register_Domain(D const & d, DS ... ds) {
    int l_pointer = sizeof...(DS);
    logic_grid_.x0[l_pointer] = d.first();
    logic_grid_.x1[l_pointer] = d.first() + d.size();
}

template <int N_RANK> 
grid_info<N_RANK> Pochoir<N_RANK>::get_phys_grid(void) {
    return phys_grid_;
}

/* Run the kernel functions stored in array of function pointers */
template <int N_RANK>
void Pochoir<N_RANK>::Run(int timestep) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_unroll(lcm_unroll_);
    timestep_ = timestep;
    /* base_case_kernel() will mimic exact the behavior of serial nested loop!
    */
    checkFlags();
    inRun = true;
    algor.base_case_kernel_guard(0 + time_shift_, timestep + time_shift_, logic_grid_, sz_pgk_, pgk_);
    inRun = false;

}

/* obase for interior and ExecSpec for boundary */
template <int N_RANK> 
void Pochoir<N_RANK>::Run_Obase(int timestep) {
    int l_total_points = 1;
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    algor.set_pgk(sz_pgk_, opgk_);
    algor.set_unroll(lcm_unroll_);
    /* this version uses 'f' to compute interior region, 
     * and 'bf' to compute boundary region
     */
    timestep_ = timestep;
    checkFlags();
    // cutting based on shorter bar
    algor.adaptive_bicut_p(0 + time_shift_, timestep + time_shift_, logic_grid_);
}


/* obase for zero-padded area! */
template <int N_RANK> template <typename F>
void Pochoir<N_RANK>::Run_Obase(int timestep, F const & f) {
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    timestep_ = timestep;
    checkFlags();
#if BICUT
#if 0
    fprintf(stderr, "Call obase_bicut\n");
    algor.obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#else
//     fprintf(stderr, "Call shorter_duo_sim_obase_bicut\n");
   // algor.sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
    algor.shorter_duo_sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
    // algor.duo_sim_obase_bicut(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#if STAT
    for (int i = 1; i < SUPPORT_RANK; ++i) {
        fprintf(stderr, "sim_count_cut[%d] = %ld\n", i, algor.sim_count_cut[i].get_value());
    }
#endif
#endif
#else
    algor.obase_m(0+time_shift_, timestep+time_shift_, logic_grid_, f);
#endif
}

/* obase for interior and ExecSpec for boundary */
template <int N_RANK> template <typename F, typename BF>
void Pochoir<N_RANK>::Run_Obase(int timestep, F const & f, BF const & bf) {
    int l_total_points = 1;
    Algorithm<N_RANK> algor(slope_);
    algor.set_phys_grid(phys_grid_);
    algor.set_thres(arr_type_size_);
    /* this version uses 'f' to compute interior region, 
     * and 'bf' to compute boundary region
     */
    timestep_ = timestep;
    checkFlags();
#if BICUT
#if 0
    fprintf(stderr, "Call obase_bicut_boundary_P\n");
    algor.obase_bicut_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#else
//    fprintf(stderr, "Call sim_obase_bicut_P\n");
//    hyper-space cut
    // algor.sim_obase_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
    // cutting based on shorter bar
    algor.shorter_duo_sim_obase_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
    // cutting based on longer bar
    // algor.duo_sim_obase_bicut_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
    // serial space cut
    // algor.obase_bicut_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#if STAT
    for (int i = 1; i < SUPPORT_RANK; ++i) {
        fprintf(stderr, "sim_count_cut[%d] = %ld\n", i, algor.sim_count_cut[i].get_value());
    }
#endif
#endif
#else
    algor.obase_boundary_p(0+time_shift_, timestep+time_shift_, logic_grid_, f, bf);
#endif
}

#endif
