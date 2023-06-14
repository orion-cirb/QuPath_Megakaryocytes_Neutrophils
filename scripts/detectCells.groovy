//------------------------------------------- PARAMETERS TO TUNE -----------------------------------------------------//
// MEGAKARYOCYTES
// Size filtering
def megakaMinArea = 400 // in µm2
def megakaMaxArea = 2000
// Intensity filtering
def megakaMinInt = 1000
def megakaMaxInt = 20000

// NEUTROPHILS
detectNeutro = true // set to true if you want to detect neutrophils, false otherwise
// Size filtering
def neutroMinArea = 20
def neutroMaxArea = 80
// Intensity filtering
def neutroMinInt = 2000
def neutroMaxInt = 8000



//----------------------------------------------- PIPELINE -----------------------------------------------------------//
import qupath.lib.gui.dialogs.Dialogs
import static qupath.lib.scripting.QP.*
import qupath.lib.objects.classes.PathClass
import qupath.opencv.ops.ImageOps
import qupath.ext.biop.cellpose.Cellpose2D
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics


// Init project
def project = getProject()

// Pixel classifier for tissue segmentation
if (!project.getPixelClassifiers().contains('Tissue')) {
    Dialogs.showErrorMessage('Problem', 'No Tissue classifier found')
    return
}
def tissueClassifier = project.getPixelClassifiers().get('Tissue')
def tissueClass = PathClass.fromString('Tissue')

// Cellpose model for megakaryocytes detection
def megakaClass = PathClass.fromString('Megakaryocyte')
megakaClass.setColor(150, 255, 0)
def megakaCellpose = buildCellposeModel('cyto2', 'EGFP', 80, true, megakaClass)

// Cellpose model for neutrophils detection
def neutroClass = PathClass.fromString('Neutrophil')
neutroClass.setColor(255,150,0)
def neutroCellpose = buildCellposeModel('cyto2', 'Cy5', 20, true, neutroClass)

// Create results file and write headers
def imageDir = new File(project.getImageList()[0].getURIs()[0]).getParent()
def resultsDir = buildFilePath(imageDir, '/Results')
if (!fileExists(resultsDir)) mkdirs(resultsDir)
def resultsFile = new File(buildFilePath(resultsDir, 'detailedResults.xls'))
resultsFile.createNewFile()
resultsFile.write("Image name\tAnnotation name\tMegakaryocyte area (um2)\tNb neutrophils\t" +
                    "Neutrophils total area (um2)\tNeutrophils mean area (um2)\tNeutrophils area std\n")
def globalResultsFile = new File(buildFilePath(resultsDir, 'globalResults.xls'))
globalResultsFile.createNewFile()
globalResultsFile.write("Image name\tAnnotation name\tAnnotation area (um2)\tTissue area (um2)\tNb megakaryocytes\t" +
        "Megakaryocytes total area (um2)\tMegakaryocytes mean area (um2)\tMegakaryocytes area std\tNb neutrophils\t" +
        "Neutrophils total area (um2)\tNeutrophils mean area (um2)\tNeutrophils area std\n")

// Loop over images in project
for (entry in project.getImageList()) {
    def imageData = entry.readImageData()
    def server = imageData.getServer()
    def cal = server.getPixelCalibration()
    def pixelWidth = cal.getPixelWidth().doubleValue()
    def imgName = entry.getImageName()
    def cy5Channel = (server.nChannels() > 2)
    setBatchProjectAndImage(project, imageData)
    setImageType('FLUORESCENCE')
    println ''
    println '------ ANALYZING IMAGE ' + imgName + ' ------'

    // Find annotations
    def annotations = getAnnotationObjects()
    if (annotations.isEmpty()) {
        Dialogs.showErrorMessage("Problem", "Please create annotations to analyze in image " + imgName)
        continue
    }
    def index = 0
    for (an in annotations) {
        if (an.getName() == null) {
            index++
            an.setName("Region_" + index)
        }
    }

    println '--- Detecting tissue with pixel classifier ---'
    deselectAll()
    selectObjects(annotations)
    createAnnotationsFromPixelClassifier(tissueClassifier, 2000, 500)
    def tissues = getAnnotationObjects().findAll({ it.getPathClass() == tissueClass })

    println '--- Detecting megakaryocytes with Cellpose ---'
    deselectAll()
    selectObjects(tissues)
    megakaCellpose.detectObjects(imageData, getSelectedObjects())
    def allMegakaryocytes = getDetectionObjects().findAll {
        it.getROI().getScaledArea(pixelWidth, pixelWidth) > megakaMinArea
                && it.getROI().getScaledArea(pixelWidth, pixelWidth) < megakaMaxArea
                && it.getMeasurementList().get('EGFP: Mean') > megakaMinInt
                && it.getMeasurementList().get('EGFP: Mean') < megakaMaxInt
    }
    allMegakaryocytes = allMegakaryocytes.collect {
        return PathObjects.createAnnotationObject(it.getROI(), it.getPathClass(), it.getMeasurementList())
    }
    clearDetections()
    addObjects(allMegakaryocytes)

    def allNeutrophils
    if (detectNeutro && cy5Channel) {
        println '--- Detecting neutrophils with Cellpose ---'
        deselectAll()
        selectObjects(tissues)
        neutroCellpose.detectObjects(imageData, getSelectedObjects())
        allNeutrophils = getDetectionObjects().findAll {
            it.getROI().getScaledArea(pixelWidth, pixelWidth) > neutroMinArea
                    && it.getROI().getScaledArea(pixelWidth, pixelWidth) < neutroMaxArea
                    && it.getMeasurementList().get('Cy5: Mean') > neutroMinInt
                    && it.getMeasurementList().get('Cy5: Mean') < neutroMaxInt
        }
        allNeutrophils = allNeutrophils.collect {
            return PathObjects.createAnnotationObject(it.getROI(), it.getPathClass(), it.getMeasurementList())
        }
        clearDetections()
        addObjects(allNeutrophils)
        resolveHierarchy()
    }

    for (an in annotations) {
        println '--- Saving results for annotation ' + an.getName() + ' ---'
        def tissue = tissues.find({an.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())})
        tissue.setName(an.getName())

        def megakaryocytes = allMegakaryocytes.findAll({tissue.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())})
        println 'Nb megakaryocytes detected = ' + megakaryocytes.size()

        def neutrophils = []
        if (detectNeutro && cy5Channel) {
            neutrophils = allNeutrophils.findAll({
                tissue.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())
                        && it.getParent() != null
                        && it.getParent().getPathClass() == megakaClass
            })
            println 'Nb neutrophils detected = ' + neutrophils.size()
        }

        deselectAll()
        selectObjects(an)
        runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons": 0.1,  "lineCap": "Round",  "removeInterior": false,  "constrainToParent": false}');
        def anDil = getAnnotationObjects().last()
        removeObject(an, true)

        deselectAll()
        selectObjects(tissue)
        runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons": 0.05,  "lineCap": "Round",  "removeInterior": false,  "constrainToParent": false}');
        def tissueDil = getAnnotationObjects().last()
        removeObject(tissue, true)

        // Save annotation, tissue and DAB cells
        tissueDil.addChildObjects(megakaryocytes)
        anDil.addChildObject(tissueDil)
        fireHierarchyUpdate()
        clearAllObjects()
        addObject(anDil)
        saveAnnotations(buildFilePath(resultsDir, imgName + "_" + anDil.getName()))

        // Write results in files
        for (megaka in megakaryocytes) {
            def neutroNb = megaka.nChildObjects()
            def neutroStats = getObjectsAreaStatistics(megaka.getChildObjects(), pixelWidth)
            resultsFile << imgName + '\t' + anDil.getName() + '\t' + megaka.getROI().getScaledArea(pixelWidth, pixelWidth) + '\t' +
                    neutroNb + '\t' + neutroStats[0] + '\t' + neutroStats[1] + '\t' + neutroStats[2] + '\n'
        }

        def megakaStats = getObjectsAreaStatistics(megakaryocytes, pixelWidth)
        def neutroStats = getObjectsAreaStatistics(neutrophils, pixelWidth)
        globalResultsFile << imgName + '\t' + anDil.getName() + '\t' + anDil.getROI().getScaledArea(pixelWidth, pixelWidth) + '\t' +
                tissueDil.getROI().getScaledArea(pixelWidth, pixelWidth) + '\t' + megakaryocytes.size() + '\t' + megakaStats[0] + '\t' +
                megakaStats[1] + '\t' + megakaStats[2] + '\t' + neutrophils.size() + '\t' + neutroStats[0] + '\t' +
                neutroStats[1] + '\t' + neutroStats[2] + '\n'
    }
    clearAllObjects()
}
println '--- All done! ---'


//------------------------------------------------- UTILS ------------------------------------------------------------//
// Build Cellpose model
def buildCellposeModel(pathModel, channel, diameter, constrainToParent, cellClass) {
    return Cellpose2D.builder(pathModel)
            .preprocess(ImageOps.Filters.median(1))     // List of preprocessing ImageOps to run on the images before exporting them
            .normalizePercentiles(1, 99)                // Percentile normalization
            .pixelSize(0.5)                     // Resolution for detection
            .channels(channel)                          // Select detection channel(s)
//          .maskThreshold(-0.2)                        // Threshold for the mask detection, defaults to 0.0
//          .flowThreshold(0.5)                         // Threshold for the flows, defaults to 0.4
            .diameter(diameter)                         // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
            .constrainToParent(constrainToParent)
            .measureShape()                             // Add shape measurements
            .measureIntensity()                         // Add intensity measurements (in all compartments)
            .classify(cellClass)
//          .doLog()
            .build()
}

// Save annotations
def saveAnnotations(imgName) {
    def path = buildFilePath(imgName + '.annot')
    def annotations = getAnnotationObjects()
    new File(path).withObjectOutputStream {
        it.writeObject(annotations)
    }
    println('Results saved')
}


// Get objects area statistics (sum, mean, std)
def getObjectsAreaStatistics(objects, pixelWidth) {
    def paramsValues = [0, 0, 0]
    if (objects.size() > 0) {
        def params = new DescriptiveStatistics()
        for (obj in objects) {
            params.addValue(obj.getROI().getScaledArea(pixelWidth, pixelWidth))
        }
        paramsValues = [params.sum, params.mean, params.standardDeviation]
    }
    return paramsValues
}
