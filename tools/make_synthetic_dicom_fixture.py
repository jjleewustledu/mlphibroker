"""Heterogeneous synthetic DICOM tree: two fabricated subjects, dynamic PET + CT + RTSTRUCT.

Every identifying value is fabricated. Layout mimics an export dump:
  root/<ExportFolder>/<mixed files>  with PET, CT and RTSTRUCT of one study in one folder, and a
  second subject whose PET and CT are split across two folders, plus a stray non-DICOM file.
"""

import datetime
import sys
from pathlib import Path

import numpy as np
import pydicom
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

root = Path(sys.argv[1])
rng = np.random.default_rng(1)


def base_ds(sop_class, patient, study, acq, series_uid, series_no, modality, desc):
    fm = FileMetaDataset()
    fm.MediaStorageSOPClassUID = sop_class
    fm.MediaStorageSOPInstanceUID = generate_uid()
    fm.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = FileDataset(None, {}, file_meta=fm, preamble=b"\0" * 128)
    ds.SOPClassUID = sop_class
    ds.SOPInstanceUID = fm.MediaStorageSOPInstanceUID
    ds.PatientName, ds.PatientID, ds.PatientBirthDate, ds.PatientSex = patient
    ds.PatientWeight = "70.5"
    ds.StudyInstanceUID, ds.FrameOfReferenceUID = study
    ds.SeriesInstanceUID = series_uid
    ds.StudyDate = acq.strftime("%Y%m%d")
    ds.SeriesDate = ds.StudyDate
    ds.AcquisitionDate = ds.StudyDate
    ds.StudyTime = acq.strftime("%H%M%S")
    ds.SeriesTime = ds.StudyTime
    ds.Modality = modality
    ds.Manufacturer = "SyntheticVendor"
    ds.ManufacturerModelName = "SynthScanner"
    ds.StationName = "SYNSTATION1"
    ds.InstitutionName = "Synthetic Institute"
    ds.DeviceSerialNumber = "SN123456789"
    ds.StudyID = "S1"
    ds.AccessionNumber = "ACC0001"
    ds.SeriesNumber = series_no
    ds.SeriesDescription = desc
    ds.ReferringPhysicianName = "Nobody^Dr"
    ds.is_little_endian = True
    ds.is_implicit_VR = False
    return ds


def write_pet(folder, patient, study, acq, nx=16, ny=16, nz=4):
    series_uid = generate_uid()
    starts, durs = [0, 30, 60, 120], [30, 30, 60, 60]
    sops = []
    for f in range(4):
        for z in range(nz):
            ds = base_ds(pydicom.uid.PositronEmissionTomographyImageStorage, patient, study, acq, series_uid, 3, "PT", "Synthetic FDG dynamic")
            t = acq + datetime.timedelta(seconds=starts[f])
            ds.AcquisitionTime = t.strftime("%H%M%S.000000")
            ds.ContentDate = ds.StudyDate
            ds.ContentTime = ds.AcquisitionTime
            ds.AcquisitionNumber = f + 1
            ds.InstanceNumber = f * nz + z + 1
            ds.ImageType = ["ORIGINAL", "PRIMARY"]
            ds.PatientPosition = "HFS"
            ds.BodyPartExamined = "BRAIN"
            ds.ImagePositionPatient = [-8.0, -8.0, float(z * 2.0)]
            ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
            ds.PixelSpacing = [1.0, 1.0]
            ds.SliceThickness = 2.0
            ds.SliceLocation = float(z * 2.0)
            ds.Rows, ds.Columns = ny, nx
            ds.SamplesPerPixel = 1
            ds.PhotometricInterpretation = "MONOCHROME2"
            ds.BitsAllocated = ds.BitsStored = 16
            ds.HighBit = 15
            ds.PixelRepresentation = 0
            ds.RescaleIntercept = 0
            ds.RescaleSlope = 1.0
            ds.Units = "BQML"
            ds.DecayCorrection = "START"
            ds.CorrectedImage = ["DECY", "ATTN", "SCAT", "RAN", "NORM"]
            ds.FrameReferenceTime = float(starts[f] * 1000 + durs[f] * 500)
            ds.ActualFrameDuration = durs[f] * 1000
            ds.ImageIndex = f * nz + z + 1
            ds.NumberOfSlices = nz
            ds.NumberOfTimeSlices = 4
            ds.ReconstructionMethod = "OSEM3D 4i21s"
            ds.AttenuationCorrectionMethod = "CT-based"
            rp = Dataset()
            rp.Radiopharmaceutical = "Fludeoxyglucose F18"
            rp.RadiopharmaceuticalStartTime = (acq - datetime.timedelta(seconds=131)).strftime("%H%M%S.000000")
            rp.RadionuclideTotalDose = 185000000.0
            rp.RadionuclideHalfLife = 6586.2
            rp.RadionuclidePositronFraction = 0.9686
            ds.RadiopharmaceuticalInformationSequence = [rp]
            ds.PixelData = (rng.integers(0, 2000, size=(ny, nx)) + f * 100).astype(np.uint16).tobytes()
            ds.save_as(folder / f"PT_{f:02d}_{z:02d}.dcm", write_like_original=False)
            if f == 0:
                sops.append((ds.SOPClassUID, ds.SOPInstanceUID, float(z * 2.0)))
    return series_uid, sops


def write_static_pet(folder, patient, study, acq, nx=16, ny=16, nz=4):
    """One-frame delayed PET of the same injection, 60 min after the dynamic scan."""
    series_uid = generate_uid()
    for z in range(nz):
        ds = base_ds(pydicom.uid.PositronEmissionTomographyImageStorage, patient, study, acq, series_uid, 4, "PT", "Synthetic FDG delayed static")
        ds.AcquisitionTime = acq.strftime("%H%M%S.000000")
        ds.ContentDate = ds.StudyDate
        ds.ContentTime = ds.AcquisitionTime
        ds.AcquisitionNumber = 1
        ds.InstanceNumber = z + 1
        ds.ImageType = ["ORIGINAL", "PRIMARY"]
        ds.PatientPosition = "HFS"
        ds.BodyPartExamined = "BRAIN"
        ds.ImagePositionPatient = [-8.0, -8.0, float(z * 2.0)]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.PixelSpacing = [1.0, 1.0]
        ds.SliceThickness = 2.0
        ds.SliceLocation = float(z * 2.0)
        ds.Rows, ds.Columns = ny, nx
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated = ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 0
        ds.RescaleIntercept = 0
        ds.RescaleSlope = 1.0
        ds.Units = "BQML"
        ds.DecayCorrection = "START"
        ds.CorrectedImage = ["DECY", "ATTN", "SCAT", "RAN", "NORM"]
        ds.FrameReferenceTime = 300.0 * 1000 / 2
        ds.ActualFrameDuration = 300 * 1000
        ds.ImageIndex = z + 1
        ds.NumberOfSlices = nz
        ds.NumberOfTimeSlices = 1
        ds.ReconstructionMethod = "OSEM3D 4i21s"
        ds.AttenuationCorrectionMethod = "CT-based"
        rp = Dataset()
        rp.Radiopharmaceutical = "Fludeoxyglucose F18"
        rp.RadiopharmaceuticalStartTime = (acq - datetime.timedelta(seconds=131 + 3600)).strftime("%H%M%S.000000")
        rp.RadionuclideTotalDose = 185000000.0
        rp.RadionuclideHalfLife = 6586.2
        rp.RadionuclidePositronFraction = 0.9686
        ds.RadiopharmaceuticalInformationSequence = [rp]
        ds.PixelData = (rng.integers(0, 500, size=(ny, nx))).astype(np.uint16).tobytes()
        ds.save_as(folder / f"PTSTATIC_{z:02d}.dcm", write_like_original=False)
    return series_uid


def write_secondary_capture(folder, patient, study, acq):
    """A dose-report style secondary capture: text would be burned into the pixels; never converted."""
    ds = base_ds(pydicom.uid.SecondaryCaptureImageStorage, patient, study, acq, generate_uid(), 502, "CT", "Patient Protocol")
    ds.ImageType = ["DERIVED", "SECONDARY", "DOSE_INFO"]
    ds.InstanceNumber = 1
    ds.Rows, ds.Columns = 32, 32
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = ds.BitsStored = 8
    ds.HighBit = 7
    ds.PixelRepresentation = 0
    ds.PixelData = (rng.integers(0, 255, size=(32, 32))).astype(np.uint8).tobytes()
    ds.save_as(folder / "SC_dose_report.dcm", write_like_original=False)


def write_ct(folder, patient, study, acq, nx=16, ny=16, nz=4):
    series_uid = generate_uid()
    for z in range(nz):
        ds = base_ds(pydicom.uid.CTImageStorage, patient, study, acq, series_uid, 2, "CT", "Synthetic CT AC")
        ds.AcquisitionTime = acq.strftime("%H%M%S.000000")
        ds.InstanceNumber = z + 1
        ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
        ds.PatientPosition = "HFS"
        ds.ImagePositionPatient = [-8.0, -8.0, float(z * 2.0)]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        ds.PixelSpacing = [1.0, 1.0]
        ds.SliceThickness = 2.0
        ds.SliceLocation = float(z * 2.0)
        ds.KVP = 120
        ds.XRayTubeCurrent = 50
        ds.ConvolutionKernel = "B30f"
        ds.Rows, ds.Columns = ny, nx
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated = ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 1
        ds.RescaleIntercept = -1024
        ds.RescaleSlope = 1.0
        ds.PixelData = rng.integers(0, 3000, size=(ny, nx)).astype(np.int16).tobytes()
        ds.save_as(folder / f"CT_{z:02d}.dcm", write_like_original=False)
    return series_uid


def write_rtstruct(folder, patient, study, acq, pet_series_uid, sops):
    ds = base_ds(pydicom.uid.RTStructureSetStorage, patient, study, acq, generate_uid(), 100, "RTSTRUCT", "MIM ROIs")
    ds.StructureSetLabel = "MIM_ROIs"
    ds.StructureSetDate = ds.StudyDate
    ds.StructureSetTime = ds.StudyTime
    ds.InstanceNumber = 1
    ref_frame = Dataset()
    ref_frame.FrameOfReferenceUID = study[1]
    rt_study = Dataset()
    rt_study.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1"
    rt_study.ReferencedSOPInstanceUID = study[0]
    rt_series = Dataset()
    rt_series.SeriesInstanceUID = pet_series_uid
    rt_series.ContourImageSequence = []
    for cls, uid, _ in sops:
        ci = Dataset()
        ci.ReferencedSOPClassUID = cls
        ci.ReferencedSOPInstanceUID = uid
        rt_series.ContourImageSequence.append(ci)
    rt_study.RTReferencedSeriesSequence = [rt_series]
    ref_frame.RTReferencedStudySequence = [rt_study]
    ds.ReferencedFrameOfReferenceSequence = [ref_frame]
    ds.StructureSetROISequence, ds.ROIContourSequence, ds.RTROIObservationsSequence = [], [], []
    rois = [(1, "LeftCaudate", [255, 0, 0], (-4.0, -4.0, 3.0)), (2, "WholeBrain", [0, 255, 0], (0.0, 0.0, 6.0))]
    for num, name, color, (cx, cy, r) in rois:
        s = Dataset()
        s.ROINumber = num
        s.ReferencedFrameOfReferenceUID = study[1]
        s.ROIName = name
        s.ROIGenerationAlgorithm = "MANUAL"
        ds.StructureSetROISequence.append(s)
        c = Dataset()
        c.ROIDisplayColor = color
        c.ReferencedROINumber = num
        c.ContourSequence = []
        for cls, uid, zpos in sops:
            cont = Dataset()
            ci = Dataset()
            ci.ReferencedSOPClassUID = cls
            ci.ReferencedSOPInstanceUID = uid
            cont.ContourImageSequence = [ci]
            cont.ContourGeometricType = "CLOSED_PLANAR"
            ang = np.linspace(0, 2 * np.pi, 13)[:-1]
            pts = []
            for a in ang:
                pts += [cx + r * np.cos(a), cy + r * np.sin(a), zpos]
            cont.NumberOfContourPoints = len(ang)
            cont.ContourData = [float(v) for v in pts]
            c.ContourSequence.append(cont)
        ds.ROIContourSequence.append(c)
        o = Dataset()
        o.ObservationNumber = num
        o.ReferencedROINumber = num
        o.ROIObservationLabel = name
        o.RTROIInterpretedType = "ORGAN"
        o.ROIInterpreter = "Nobody^Dr"
        ds.RTROIObservationsSequence.append(o)
    ds.save_as(folder / "RS_MIM.dcm", write_like_original=False)


# Subject A: everything in one export folder
a_dir = root / "export_batch1" / "A_mixed"
a_dir.mkdir(parents=True, exist_ok=True)
pA = ("Synthetic^Alpha", "SYN000001", "19700101", "O")
sA = (generate_uid(), generate_uid())
acqA = datetime.datetime(2001, 2, 3, 12, 0, 11)
pet_uid, sops = write_pet(a_dir, pA, sA, acqA)
write_static_pet(a_dir, pA, sA, acqA + datetime.timedelta(minutes=60))
write_secondary_capture(a_dir, pA, sA, acqA)
write_ct(a_dir, pA, sA, acqA - datetime.timedelta(minutes=5))
write_rtstruct(a_dir, pA, sA, acqA, pet_uid, sops)
(a_dir / "notes.txt").write_text("stray non-DICOM file\n")
# Subject B: PET and CT in separate folders, two studies (two visits)
pB = ("Synthetic^Beta", "SYN000002", "19800202", "O")
for k, day in enumerate((5, 19)):
    sB = (generate_uid(), generate_uid())
    acqB = datetime.datetime(2001, 3, day, 9, 30, 0)
    (root / "export_batch2" / f"B_visit{k}" / "PET").mkdir(parents=True, exist_ok=True)
    (root / "export_batch2" / f"B_visit{k}" / "CT").mkdir(parents=True, exist_ok=True)
    write_pet(root / "export_batch2" / f"B_visit{k}" / "PET", pB, sB, acqB)
    write_ct(root / "export_batch2" / f"B_visit{k}" / "CT", pB, sB, acqB)
print("heterogeneous fixture written under", root)
