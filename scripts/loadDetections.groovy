import org.apache.commons.io.FilenameUtils
import static qupath.lib.scripting.QP.*


// Get image directory
def imageDir = new File(project.getImageList()[0].getUris()[0]).getParent()
def imageName = getCurrentImageData().getServer().getMetadata().getName()

// Delete all annotations
clearAllObjects()

// Load annotations files for current image
def p = ~/${imageName}.*\.annot/
def resultsDir = new File(buildFilePath(imageDir+'/Results'))
resultsDir.eachFileMatch(p) {file ->
    new File(file.path).withObjectInputStream {
        def annotations = it.readObject()
        print('Adding annotation ' + annotations.toString())
        addObjects(annotations)
    }
}
resolveHierarchy()