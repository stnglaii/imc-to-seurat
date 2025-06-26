slide.seq = CreateSeuratObject(counts = COUNTS_MTX, assay="Spatial")

coords <- FetchData(seurat_object, vars = c("image_roi", "cell_pos_x", "cell_pos_y")) %>% 
  filter(image_roi == 1) %>% 
  select(-image_roi)


seurat_object@images$image =  new(
  Class = 'SlideSeq',
  assay = "RNA",
  key = "image_",
  coordinates = coords
)

seurat_object@images[["ROI_1"]] <- CreateFOV(coords, type = "centroids", assay = "RNA")




Read10X_Image <- function(
    image.dir,
    image.name = "tissue_lowres_image.png",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE,
    image.type = "VisiumV2"
) {
  # Validate the `image.type` parameter.
  image.type <- match.arg(image.type, choices = c("VisiumV1", "VisiumV2"))
  
  # Read in the H&E stain image.
  image <- png::readPNG(
    source = file.path(
      image.dir,
      image.name
    )
  )
  
  # Read in the scale factors.
  scale.factors <- Read10X_ScaleFactors(
    filename = file.path(image.dir, "scalefactors_json.json")
  )
  
  # Read in the tissue coordinates as a data.frame.
  coordinates <- Read10X_Coordinates(
    filename = Sys.glob(file.path(image.dir, "*tissue_positions*")),
    filter.matrix
  )
  
  # Use the `slice` value to populate a Seurat-style identifier for the image.
  key <- Key(slice, quiet = TRUE)
  
  # Return the specified `image.type`.
  if (image.type == "VisiumV1") {
    visium.v1 <- new(
      Class = image.type,
      assay = assay,
      key = key,
      coordinates = coordinates,
      scale.factors = scale.factors,
      image = image
    )
    
    # As of v5.1.0 `Radius.VisiumV1` no longer returns the value of the 
    # `spot.radius` slot and instead calculates the value on the fly, but we 
    # can populate the static slot in case it's depended on.
    visium.v1@spot.radius <- Radius(visium.v1)
    
    return(visium.v1)
  }
  
  # If `image.type` is not "VisiumV1" then it must be "VisiumV2".
  stopifnot(image.type == "VisiumV2")
  
  # Create an `sp` compatible `FOV` instance.
  fov <- CreateFOV(
    coordinates[, c("imagerow", "imagecol")],
    type = "centroids",
    radius = scale.factors[["spot"]],
    assay = assay,
    key = key
  )
  
  # Build the final `VisiumV2` instance, essentially just adding `image` and
  # `scale.factors` to the `fov`.
  visium.v2 <- new(
    Class = "VisiumV2",
    boundaries = fov@boundaries,
    molecules = fov@molecules,
    assay = fov@assay,
    key = fov@key,
    image = image,
    scale.factors = scale.factors
  )
  
  return(visium.v2)
}