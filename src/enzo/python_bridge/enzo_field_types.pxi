cdef extern from "typedefs.h" nogil:

    enum: Density
    enum: TotalEnergy
    enum: InternalEnergy
    enum: Pressure
    enum: Velocity1
    enum: Velocity2
    enum: Velocity3
    enum: ElectronDensity
    enum: HIDensity
    enum: HIIDensity
    enum: HeIDensity
    enum: HeIIDensity
    enum: HeIIIDensity
    enum: HMDensity
    enum: H2IDensity
    enum: H2IIDensity
    enum: DIDensity
    enum: DIIDensity
    enum: HDIDensity
    enum: SNColour
    enum: Metallicity
    enum: ExtraType0
    enum: ExtraType1
    enum: kphHI
    enum: PhotoGamma
    enum: kphHeI
    enum: gammaHeI
    enum: kphHeII
    enum: gammaHeII
    enum: kdissH2I
    enum: GravPotential
    enum: Acceleration0
    enum: Acceleration1
    enum: Acceleration2
    enum: RadPressure0
    enum: RadPressure1
    enum: RadPressure2
    enum: Emissivity0
  
    enum: gParticlePosition
    enum: gParticleVelocity
    enum: gParticleMass
    enum: gParticleAcceleration
    enum: gParticleNumber
    enum: gParticleType
    enum: gParticleAttribute
    enum: gPotentialField
    enum: gAccelerationField
    enum: gGravitatingMassField
    enum: gFlaggingField
    enum: gVelocity
  
    enum: Bfield1
    enum: Bfield2
    enum: Bfield3
    enum: PhiField
    enum: Phi_pField
    enum: DebugField
  
    enum: DrivingField1
    enum: DrivingField2
    enum: DrivingField3
  
    enum: AccelerationField1
    enum: AccelerationField2
    enum: AccelerationField3
  
    enum: Galaxy1Colour
    enum: Galaxy2Colour
  
    enum: Mach
    enum: PreShockTemperature
    enum: PreShockDensity
    enum: CRDensity
  
    enum: CIDensity
    enum: CIIDensity
    enum: OIDensity
    enum: OIIDensity
    enum: SiIDensity
    enum: SiIIDensity
    enum: SiIIIDensity
    enum: CHIDensity
    enum: CH2IDensity
    enum: CH3IIDensity
    enum: C2IDensity
    enum: COIDensity
    enum: HCOIIDensity
    enum: OHIDensity
    enum: H2OIDensity
    enum: O2IDensity
  
    enum: MBHColour
    enum: ForbiddenRefinement
  
  
    enum: RadiationFreq0
    enum: RadiationFreq1
    enum: RadiationFreq2
    enum: RadiationFreq3
    enum: RadiationFreq4
    enum: RadiationFreq5
    enum: RadiationFreq6
    enum: RadiationFreq7
    enum: RadiationFreq8
    enum: RadiationFreq9
  
    enum: RaySegments
  
    enum: FieldUndefined

cdef field_enums = {
    "Density" : Density,
    "TotalEnergy" : TotalEnergy,
    "GasEnergy" : InternalEnergy,
    "Pressure" : Pressure,
    "x-velocity" : Velocity1,
    "y-velocity" : Velocity2,
    "z-velocity" : Velocity3,
    "Electron_Density" : ElectronDensity,
    "HI_Density" : HIDensity,
    "HII_Density" : HIIDensity,
    "HeI_Density" : HeIDensity,
    "HeII_Density" : HeIIDensity,
    "HeIII_Density" : HeIIIDensity,
    "HM_Density" : HMDensity,
    "H2I_Density" : H2IDensity,
    "H2II_Density" : H2IIDensity,
    "DI_Density" : DIDensity,
    "DII_Density" : DIIDensity,
    "HDI_Density" : HDIDensity,
    "SNColour" : SNColour,
    "Metallicity" : Metallicity,
    "ExtraType0" : ExtraType0,
    "ExtraType1" : ExtraType1,
    "kphHI" : kphHI,
    "PhotoGamma" : PhotoGamma,
    "kphHeI" : kphHeI,
    "gammaHeI" : gammaHeI,
    "kphHeII" : kphHeII,
    "gammaHeII" : gammaHeII,
    "kdissH2I" : kdissH2I,
    "GravPotential" : GravPotential,
    "Acceleration0" : Acceleration0,
    "Acceleration1" : Acceleration1,
    "Acceleration2" : Acceleration2,
    "RadPressure0" : RadPressure0,
    "RadPressure1" : RadPressure1,
    "RadPressure2" : RadPressure2,
    "Emissivity0" : Emissivity0,
  
    "gParticlePosition" : gParticlePosition,
    "gParticleVelocity" : gParticleVelocity,
    "gParticleMass" : gParticleMass,
    "gParticleAcceleration" : gParticleAcceleration,
    "gParticleNumber" : gParticleNumber,
    "gParticleType" : gParticleType,
    "gParticleAttribute" : gParticleAttribute,
    "gPotentialField" : gPotentialField,
    "gAccelerationField" : gAccelerationField,
    "gGravitatingMassField" : gGravitatingMassField,
    "gFlaggingField" : gFlaggingField,
    "gVelocity" : gVelocity,
  
    "Bfield1" : Bfield1,
    "Bfield2" : Bfield2,
    "Bfield3" : Bfield3,
    "PhiField" : PhiField,
    "Phi_pField" : Phi_pField,
    "DebugField" : DebugField,
  
    "DrivingField1" : DrivingField1,
    "DrivingField2" : DrivingField2,
    "DrivingField3" : DrivingField3,
  
    "AccelerationField1" : AccelerationField1,
    "AccelerationField2" : AccelerationField2,
    "AccelerationField3" : AccelerationField3,
  
    "Galaxy1Colour" : Galaxy1Colour,
    "Galaxy2Colour" : Galaxy2Colour,
  
    "Mach" : Mach,
    "PreShockTemperature" : PreShockTemperature,
    "PreShockDensity" : PreShockDensity,
    "CRDensity" : CRDensity,
  
    "CI_Density" : CIDensity,
    "CII_Density" : CIIDensity,
    "OI_Density" : OIDensity,
    "OII_Density" : OIIDensity,
    "SiI_Density" : SiIDensity,
    "SiII_Density" : SiIIDensity,
    "SiIII_Density" : SiIIIDensity,
    "CHI_Density" : CHIDensity,
    "CH2I_Density" : CH2IDensity,
    "CH3II_Density" : CH3IIDensity,
    "C2I_Density" : C2IDensity,
    "COI_Density" : COIDensity,
    "HCOII_Density" : HCOIIDensity,
    "OHI_Density" : OHIDensity,
    "H2OI_Density" : H2OIDensity,
    "O2I_Density" : O2IDensity,
  
    "MBHColour" : MBHColour,
    "ForbiddenRefinement" : ForbiddenRefinement,
  
  
    "RadiationFreq0" : RadiationFreq0,
    "RadiationFreq1" : RadiationFreq1,
    "RadiationFreq2" : RadiationFreq2,
    "RadiationFreq3" : RadiationFreq3,
    "RadiationFreq4" : RadiationFreq4,
    "RadiationFreq5" : RadiationFreq5,
    "RadiationFreq6" : RadiationFreq6,
    "RadiationFreq7" : RadiationFreq7,
    "RadiationFreq8" : RadiationFreq8,
    "RadiationFreq9" : RadiationFreq9,
  
    "RaySegments" : RaySegments,
  
    "FieldUndefined" : FieldUndefined
}
