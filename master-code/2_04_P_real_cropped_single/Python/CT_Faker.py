import pydicom
import os
import re

def fake_pixel_data(slice_index, directory="../DICOM/", output_directory="../DICOM_autofake/"):
    """
    Kopiert die PixelData von der ausgewählten DICOM-Datei auf ihre benachbarten Dateien.
    
    Args:
        slice_index (int): Index der Datei, die gefälscht wird, z.B. 46 für IMG0046_cropped.dcm.
        directory (str): Verzeichnis der Original-DICOM-Dateien.
        output_directory (str): Verzeichnis, in das die modifizierten DICOM-Dateien gespeichert werden.
    """
    # Erstellen des Dateipfads der ausgewählten DICOM-Datei basierend auf dem Index
    fake_file = f"{directory}IMG00{slice_index}_cropped.dcm"
    
    # Prüfen, ob die ausgewählte Fake-Datei existiert
    if not os.path.exists(fake_file):
        raise FileNotFoundError(f"{fake_file} existiert nicht.")
    
    # Laden der Fake-DICOM-Datei
    fake_dataset = pydicom.dcmread(fake_file)

    # Schleife über den Bereich der Originaldateien um die ausgewählte Datei und ihre Nachbarn
    for i in range(slice_index - 1, slice_index + 2):
        original_file = f"{directory}IMG00{i}_cropped.dcm"
        
        # Prüfen, ob die Datei existiert
        if not os.path.exists(original_file):
            print(f"{original_file} existiert nicht. Überspringen...")
            continue
        
        # Laden der Original-DICOM-Datei
        original_dataset = pydicom.dcmread(original_file)
        
        # Ersetzen der PixelData im originalen Dataset mit der aus dem Fake-Dataset
        original_dataset.PixelData = fake_dataset.PixelData
        
        # Festlegen des Ausgabedateinamens: 
        # Wenn es die ausgewählte Fake-Datei ist, den ursprünglichen Namen behalten
        # Sonst den Namen mit "_fake" speichern
        if original_file == fake_file:
            output_file = f"{output_directory}IMG00{i}_cropped.dcm"
        else:
            output_file = f"{output_directory}IMG00{i}_cropped_fake.dcm"
        
        # Speichern der modifizierten Datei
        original_dataset.save_as(output_file)
        
        print(f"PixelData erfolgreich von {fake_file} nach {original_file} kopiert und gespeichert als {output_file}.")
    
    print("Die ausgewählte Datei und ihre nächsten Nachbarn wurden erfolgreich verarbeitet.")

# Beispielaufruf der Funktion mit Index 46
fake_pixel_data(46)
