# I got these boxes from Klein 2017
# who got them from https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2015GL065627
# who got them from https://link.springer.com/article/10.1007/s00382-013-1945-z

lowcloudboxes = [ # lower-left cornerof each box
    (-110, -30), # Peru
    (-155, 15), # california
    (-28, -25), # Namibia
    (-55, 10), # Canary
    (75, -38), # Australia
]
lowcloudnames = ["Peru", "California", "Namibia", "Canary", "Australia"]
low_cloud_box_lat = 20
low_cloud_box_lon = 40

low_cloud_boxes_corners = Dict(name => box for (box, name) in zip(lowcloudboxes, lowcloudnames))

low_cloud_boxes_ERA5 = Dict()

for (box, name) in zip(lowcloudboxes, lowcloudnames)
    # the area must be specified as north, west, south, east
    # for use when downloading ERA5 data from CDSAPI
    low_cloud_boxes_ERA5[name] = [box[2]+low_cloud_box_lat, box[1], box[2], box[1] + low_cloud_box_lon]
end
