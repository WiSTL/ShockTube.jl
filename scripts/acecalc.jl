### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 0bf7ad90-e88f-11ea-147b-eb2f2eb606be
using PlutoUI, ShockTube, Unitful

# ╔═╡ 2f4c7d80-e893-11ea-3771-d3c8088a873a
md"## Calculating from flow rate settings"

# ╔═╡ 1f7b19b0-e88f-11ea-27b0-8db91c84aedb
Ace = Species("Acetone");

# ╔═╡ 27f9c870-e88f-11ea-0808-4ba3e41fc037
He = Species("Helium");

# ╔═╡ a6076a9e-e890-11ea-232c-61bbfd839410
begin
	Ace_rate_slider = @bind Ace_rate NumberField(LinRange(0, 20, 1001), default=10)
	He_rate_slider = @bind He_rate NumberField(LinRange(0, 50, 1001), default=30)

	md"""
	Acetone flow rate (sccm) $Ace_rate_slider
	
	Helium flow rate (slpm)  $He_rate_slider
	"""
end

# ╔═╡ 8cc5a110-e8b8-11ea-2cec-7f83d4d4f7ba
v̇Ace = Ace_rate*u"cm^3/minute"

# ╔═╡ 96e359ce-e8b8-11ea-29ae-d3c78b2e1dba
v̇He = He_rate*u"L/minute"

# ╔═╡ 6db8dae0-e88f-11ea-1355-bbce79b13063
ṅAce = ρm(Ace) * v̇Ace |> u"mol/s"

# ╔═╡ 70091700-e890-11ea-0307-31689d3e8e9f
ṅHe = ρm(He) * v̇He |> u"mol/s"

# ╔═╡ 85bf68b0-e890-11ea-04d6-e1cfff79278d
χAce = ṅAce / (ṅAce + ṅHe)

# ╔═╡ 3e9944d0-e893-11ea-30de-332f76096802
md"## Calculating from desired flow rate & acetone concentration"

# ╔═╡ b53e58a0-e893-11ea-111e-37b14cf3dcab
begin
	v̇_slider = @bind v̇ NumberField(LinRange(0, 50, 101))
	χAce_slider = @bind χAce_in NumberField(LinRange(0, 0.25, 101)) 
	md"""
	Net volume flow rate (slpm) $v̇_slider
	
	Acetone fraction $χAce_slider
	"""
end

# ╔═╡ 692684b0-e893-11ea-2af9-756ab48bcdd0
HeAce = Mixture(["Helium", "Acetone"], zs=[1 - χAce_in, χAce_in])

# ╔═╡ 3bcce8f0-e894-11ea-0ff0-25d9c9785a15
ṅ = ρm(HeAce) * v̇*u"L/minute" |> u"mol/s"

# ╔═╡ ab9018a0-e895-11ea-10ab-03e7704de77d
Ace_slpm = ṅ * χAce_in / ρm(Ace) |> u"cm^3 / minute"

# ╔═╡ 6b7a74e0-e895-11ea-015c-aba28e723af1
He_slpm = ṅ * (1 - χAce_in) / ρm(He) |> u"L/minute"

# ╔═╡ 209a1e20-eef9-11ea-2e0a-b3b5b43861b6
md""" ## Shock jump conditions for helium/acetone mixture"""

# ╔═╡ c238a120-eef9-11ea-35cd-f744665ba19a
begin
	Mach_slider = @bind Mach NumberField(range(1.1, 2.8, step=0.01))
	md"""
	Mach number $Mach_slider
	"""
end

# ╔═╡ 1da4d2f0-eef9-11ea-0752-cf0358151391
driver_gas = Species("N2", T=800)

# ╔═╡ 13cbc450-eef9-11ea-24e0-0fa6efd42922
driver_pressure = p(shockcalc(driver_gas, HeAce, Mach).driver) |> u"psi"

# ╔═╡ Cell order:
# ╟─2f4c7d80-e893-11ea-3771-d3c8088a873a
# ╠═0bf7ad90-e88f-11ea-147b-eb2f2eb606be
# ╠═1f7b19b0-e88f-11ea-27b0-8db91c84aedb
# ╠═27f9c870-e88f-11ea-0808-4ba3e41fc037
# ╟─a6076a9e-e890-11ea-232c-61bbfd839410
# ╟─8cc5a110-e8b8-11ea-2cec-7f83d4d4f7ba
# ╟─96e359ce-e8b8-11ea-29ae-d3c78b2e1dba
# ╟─6db8dae0-e88f-11ea-1355-bbce79b13063
# ╟─70091700-e890-11ea-0307-31689d3e8e9f
# ╟─85bf68b0-e890-11ea-04d6-e1cfff79278d
# ╟─3e9944d0-e893-11ea-30de-332f76096802
# ╟─b53e58a0-e893-11ea-111e-37b14cf3dcab
# ╠═692684b0-e893-11ea-2af9-756ab48bcdd0
# ╟─3bcce8f0-e894-11ea-0ff0-25d9c9785a15
# ╟─ab9018a0-e895-11ea-10ab-03e7704de77d
# ╟─6b7a74e0-e895-11ea-015c-aba28e723af1
# ╟─209a1e20-eef9-11ea-2e0a-b3b5b43861b6
# ╠═c238a120-eef9-11ea-35cd-f744665ba19a
# ╠═1da4d2f0-eef9-11ea-0752-cf0358151391
# ╠═13cbc450-eef9-11ea-24e0-0fa6efd42922
