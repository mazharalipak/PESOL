function [ tnr ] = tnr_init( mpc )
%init_tnr Initialize the TNR algorithm data structure
%   Detailed explanation goes here

    tnr = struct();
    
    [baseMVA, bus, gen, branch] = deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch);
    
    [tnr.ref, tnr.pv, tnr.pq] = bustypes(bus, gen);
    [tnr.Ybus, tnr.Yf, tnr.Yt] = makeYbus(baseMVA, bus, branch);
    tnr.S0 = makeSbus(baseMVA, bus, gen); % ignoring the ZIP loads for a moment
    tnr.dS = tnr.S0;
    
    tnr.npv = length(tnr.pv);
    tnr.npq = length(tnr.pq);

    % Following the runpf procedures
    
    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
        VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
        MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
        QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
    
    on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
    gbus = gen(on, GEN_BUS);                %% what buses are they at?

    V0  = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
    vcb = ones(size(V0));           %% create mask of voltage-controlled buses
    vcb(tnr.pq) = 0;                %% exclude PQ buses
    k = find(vcb(gbus));            %% in-service gens at v-c buses
    V0(gbus(k)) = gen(on(k), VG) ./ abs(V0(gbus(k))).* V0(gbus(k));
    
    tnr.V0 = V0;
    
    % Default functions, change if needed
    % Regular power flow
    tnr.F = @F_pf;
    tnr.J = @J_pf;
    tnr.dJdxl = @dJdxl_pf;
    
    % Use svd transversality condition
    tnr.G = @G_svd;
    tnr.dGdz = @dGdz_svd;
end

