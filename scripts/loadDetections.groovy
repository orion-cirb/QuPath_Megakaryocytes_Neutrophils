import org.apache.commons.io.FilenameUtils
import static qupath.lib.scripting.QP.*
import qupath.lib.gui.dialogs.Dialogs


// Get image directory
def imageName = getCurrentImageData().getServer().getMetadata().getName()
def imageDir = new File(project.getImageList()[0].getURIs()[0]).getParent()
def resultsDir = Dialogs.promptForDirectory('Select results directory', new File(imageDir))

// Delete all annotations
clearAllObjects()

// Load annotations files for current image
def p = ~/${imageName}.*\.annot/
resultsDir.eachFileMatch(p) {file ->
    new File(file.path).withObjectInputStream {
        def annotations = it.readObject()
        print('Adding annotation ' + annotations.toString())
        addObjects(annotations)
    }
}
resolveHierarchy()