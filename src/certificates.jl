export Putinar, relax!

abstract type AbstractCertificate end

mutable struct Putinar <: AbstractCertificate
    halfdegree::Int
    chordal::Bool
end

function Putinar(halfdegree; chordal=false)
    return Putinar(halfdegree, chordal)
end

function Base.show(io::IO, cert::Putinar)
    println(io, "Putinar Certificate of total degree $(2*cert.halfdegree).")
    if cert.chordal
        println(io, "Chordal sparsity used.")
    end
end

            
function relax!(m::PolyModel, cert::Putinar, solver_factory)
    model = SOSModel()

    if objective_sense(m) == :Max
        @variable model t
        @objective model Min t
        sosp = t - objective_function(m)
    elseif objective_sense(m) == :Min
        @variable model t
        @objective model Max t
        sosp = objective_function(m) - t
    else
        sosp = objective_function(m)
    end
    if cert.chordal
        chordal_putinar(sosp, 2*cert.halfdegree, constraints(m), model = model)
    else
        putinar(sosp, 2*cert.halfdegree, constraints(m), model = model)
    end
    

    time = @elapsed optimize!(model, solver_factory)
    m.relaxation[:model] = model
    m.relaxation[:cert] = cert
    m.relaxation[:solver] = solver_factory
    m.relaxation[:status] = termination_status(model)
    m.relaxation[:PSTAT] = primal_status(model)
    m.relaxation[:DSTAT] = dual_status(model)
    m.relaxation[:bound] = objective_value(model)
    m.relaxation[:time_solve] = time
    return m
end


